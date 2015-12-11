"""
Simulate images using galsim instead of ngmix
"""
from __future__ import print_function
import os
import pprint

import numpy

import ngmix
from . import sim as ngmixsim
from .sim import srandu

from .util import TryAgainError, load_gmixnd

try:
    import galsim
    have_galsim=True
except ImportError:
    have_galsim=False

class SimGS(dict):
    def __init__(self, sim_conf):
        """
        Simulate images

        example
        -------
        sim=NGMixSim(conf)

        for i in xrange(1000):
            imdict = si.get_image()
            # process image
        """

        self.update(sim_conf)

        self._setup()

        pprint.pprint(self)
        if 'seed' in self:
            print("    using seed:",self['seed'])
            numpy.random.seed(self['seed'])


    def get_image(self):
        """
        get a randomized galaxy image
        """

        gal, gal_pars, psf = self._get_galsim_objects()

        psf_obs = self._make_obs(psf, s2n=self['psf']['s2n'])
        gal_obs = self._make_obs(gal)

        s2n = self._get_expected_s2n(gal_obs.image_nonoise)
        print("    s2n expected:",s2n)

        gal_obs.set_psf(psf_obs)

        save_pars=[gal_pars['r50'], gal_pars['flux']]
        return {'obs':gal_obs,
                's2n':s2n,
                'model':self['model'],
                'gal_info':gal_pars,
                'pars':save_pars,
                'psf_obj': psf,
                'gal_obj': gal}


    def _make_obs(self, gs_obj, s2n=None):
        """
        get an ngmix Observation

        for psfs, send s2n=
        """


        nrows,ncols=self['stamp_size']
        gsimage = gs_obj.drawImage(nx=ncols,
                                   ny=nrows,
                                   scale=1.0,
                                   dtype=numpy.float64)

        im0 = gsimage.array
        if s2n is not None:
            image_nonoise, image, flux = self._scale_and_add_noise(im0, s2n)
        else:
            image_nonoise = im0.copy()
            image = self._add_noise(image_nonoise)

        jacob = self._get_jacobian(image_nonoise)

        weight = numpy.zeros( image.shape ) + self['ivar']

        obs = ngmix.Observation(image, weight=weight, jacobian=jacob)

        # monkey patching
        obs.image_nonoise = image_nonoise

        if False and s2n is None:
            self._compare_images(image_nonoise,image,label1='im',label2='noisy')

        return obs


    def _compare_images(self, im1, im2, **keys):
        import images
        keys['width']=1000
        keys['height']=1000
        images.compare_images(im1, im2, **keys)
        key=raw_input('hit a key: ')
        if key=='q':
            stop


    def _add_noise(self, im0):
        nim = numpy.random.normal(loc=0.0,
                                  scale=self['noise'],
                                  size=im0.shape)
        noisy_image = im0 + nim
        return noisy_image


    def _scale_and_add_noise(self, im0, s2n):
        """
        find the flux that gives the requested s2n.
        """

        flux = self._get_flux_from_s2n(im0, s2n)

        imsum = im0.sum()

        factor = flux/imsum

        scaled_image = im0 * factor

        nim = numpy.random.normal(loc=0.0,
                                  scale=self['noise'],
                                  size=im0.shape)
        noisy_image = scaled_image + nim

        return scaled_image, noisy_image, flux

    def _get_jacobian(self, image):
        """
        find the best center and set the jacobian center there
        """
        fitter = quick_fit_gauss(image)
        row,col = fitter.get_gmix().get_cen()
        #print("    row,col:",row,col)

        scale=1.0
        return ngmix.Jacobian(row, col,
                              scale,
                              0.0,
                              0.0,
                              scale)

    def _get_galsim_objects(self):
        """
        get the galaxy and psf galsim objects
        """

        gal0, gal_pars = self._get_gal_obj()
        psf  = self._get_psf_obj(gal_pars['cenoff'])

        gal = galsim.Convolve([psf, gal0])

        return gal, gal_pars, psf

    def _get_gal_obj(self):
        """
        get the galsim object for the galaxy model
        """

        pars=self._get_galaxy_pars()

        flux = pars['flux']
        r50 = pars['r50']

        g1,g2=pars['g']
        shear,shindex = self.shear_pdf.get_shear()
        pars['shear'] = shear
        pars['shear_index'] = shindex

        cenoff = pars['cenoff']

        if pars['model']=='gauss':
            gal = galsim.Gaussian(flux=flux, half_light_radius=r50)
        elif pars['model']=='exp':
            gal = galsim.Exponential(flux=flux, half_light_radius=r50)
        elif pars['model']=='dev':
            gal = galsim.DeVaucouleurs(flux=flux, half_light_radius=r50)
        else:
            raise ValueError("bad galaxy model: '%s'" % self['model'])

        # first give it an intrinsic shape
        gal = gal.shear(g1=g1, g2=g2)

        # now shear it
        gal = gal.shear(g1=shear.g1, g2=shear.g2)

        # in the demos, the shift was always applied after the shear, not sure
        # if it matters
        if cenoff is not None:
            gal = gal.shift(dx=cenoff[0], dy=cenoff[1])

        tup=(r50,cenoff)
        print("    r50: %g cenoff: %s" % tup)

        return gal, pars


    def _get_psf_obj(self, cenoff):
        """
        get the galsim object for the psf
        """

        pspec = self['psf']

        model=pspec['model']

        r50 = self._get_psf_r50()

        if model=='moffat':
            psf = galsim.Moffat(beta=pspec['beta'], half_light_radius=r50)
        elif model=='gauss':
            psf = galsim.Gaussian(half_light_radius=r50)
        else:
            raise ValueError("bad psf model: '%s'" % model)

        psf_g1, psf_g2 = self._get_psf_shape()
        psf = psf.shear(g1=psf_g1, g2=psf_g2)

        if cenoff is not None:
            psf = psf.shift(dx=cenoff[0], dy=cenoff[1])

        return psf

    def _get_psf_r50(self):
        if self.psf_r50_pdf is not None:
            r50 = self.psf_r50_pdf.sample()
            print("    psf (pdf) r50: %g" % r50)
        else:
            r50 = self['psf']['r50']
        return r50

    def _get_psf_shape(self):
        pspec = self['psf']

        if self.psf_ellip_pdf is not None:
            psf_g1, psf_g2 = self.psf_ellip_pdf.sample()
            print("    psf (pdf) shape: %g %g" % (psf_g1, psf_g2))

        elif 'randomized_orientation' in pspec:
            ro=pspec['randomized_orientation']
            if ro["dist"]=="uniform":
                angle = numpy.random.random()*2*numpy.pi
                psf_shape = ngmix.Shape(ro['magnitude'], 0.0)
                psf_shape.rotate(angle)
                psf_g1 = psf_shape.g1
                psf_g2 = psf_shape.g2
                print("    psf rand orient. shape: %g %g" % (psf_g1, psf_g2))
            else:
                raise ValueError("only uniform randomized psf orientation for now")
        else:
            psf_g1=pspec['shape'][0]
            psf_g2=pspec['shape'][1]

        return psf_g1, psf_g2

    def _get_galaxy_pars(self):
        """
        Get pair parameters

        if not random, then the mean pars are used, except for cen and g1,g2
        which are zero
        """

        if self.cen_pdf is not None:
            coff1 = self.cen_pdf.sample()
            coff2 = self.cen_pdf.sample()

            cenoff=(coff1,coff2)
        else:
            cenoff=None

        g1,g2 = self.g_pdf.sample2d()

        flux = self.flux_pdf.sample()
        if self.flux_is_in_log:
            flux = numpy.exp(flux)

        # this is the round r50
        r50 = self.r50_pdf.sample()

        pars = {'model':self['model'],
                'g':(g1,g2),
                'flux':flux,
                'r50':r50,
                'cenoff':cenoff}

        return pars

    '''
    def _get_galaxy_pars(self):
        """
        Get pair parameters

        if not random, then the mean pars are used, except for cen and g1,g2
        which are zero
        """

        if self.cen_pdf is not None:
            coff1 = self.cen_pdf.sample()
            coff2 = self.cen_pdf.sample()

            cenoff=(coff1,coff2)
        else:
            cenoff=None

        g1,g2 = self.g_pdf.sample2d()

        r50 = self.r50_pdf.sample()

        pars = {'model':self['model'],
                'g':(g1,g2),
                'r50':r50,
                'cenoff':cenoff}

        return pars
    '''
    def _get_flux_from_s2n(self, image, s2n):
        """
        get the flux required to give the image the requested s/n
        """
        s2n_raw = self._get_expected_s2n(image)

        flux2use = s2n/s2n_raw

        if False:
            imsum=image.sum()
            print("    image sum:",imsum)
            timage = image * flux2use/imsum
            s2n_got = self._get_expected_s2n(timage)
            print("    requested s2n: %g got %g" % (s2n, s2n_got) )
        return flux2use

    def _get_expected_s2n(self, image):
        s2n = numpy.sqrt( (image**2).sum() * self['ivar'] )
        return s2n


    def _setup(self):

        self['ivar'] = 1.0/self['noise']**2
        self['model'] = self['obj_model']['model']
        self._set_pdfs()

    def _set_pdfs(self):
        """
        Set all the priors
        """

        self._set_psf_pdf()
        self._set_flux_pdf()
        self._set_g_pdf()
        self._set_r50_pdf()
        self._set_cen_pdf()
        self._set_shear_pdf()

    def _set_shear_pdf(self):
        from .shearpdf import ConstShearGenerator
        shconf = self['shear']
        if shconf['type'] == 'const':
            pdf = ConstShearGenerator(shconf['shears'])
        else:
            raise ValueError("only shear 'const' for now")

        self.shear_pdf=pdf

    def _set_psf_pdf(self):
        pspec = self['psf']

        if isinstance(pspec['r50'], dict):
            r50pdf = pspec['r50']
            assert r50pdf['type']=="lognormal","r50 pdf log normal for now"

            self.psf_r50_pdf = ngmix.priors.LogNormal(r50pdf['mean'],
                                                      r50pdf['sigma'])
        else:
            self.psf_r50_pdf=None

        if isinstance(pspec['shape'],dict):
            ppdf=pspec['shape']
            assert ppdf['type']=="normal2d"
            self.psf_ellip_pdf=ngmix.priors.SimpleGauss2D(ppdf['cen'][0],
                                                          ppdf['cen'][1],
                                                          ppdf['sigma'][0],
                                                          ppdf['sigma'][1])
        else:
            self.psf_ellip_pdf=None

    def _set_g_pdf(self):
        g_spec=self['obj_model']['g']
        self.g_pdf=ngmix.priors.GPriorBA(g_spec['sigma'])

    def _set_r50_pdf(self):
        r50spec = self['obj_model']['r50']
        
        if r50spec['type']=='uniform':
            r50_r = r50spec['range']
            self.r50_pdf=ngmix.priors.FlatPrior(r50_r[0], r50_r[1])
        elif r50spec['type']=='lognormal':
            self.r50_pdf=ngmix.priors.LogNormal(r50spec['mean'],
                                                r50spec['sigma'])
        else:
            raise ValueError("bad r50 pdf type: '%s'" % r50spec['type'])

    def _set_cen_pdf(self):

        cr=self['obj_model']['cen_shift']
        if cr is None:
            self.cen_pdf=None
        else:
            self.cen_pdf=ngmix.priors.FlatPrior(-cr['radius'], cr['radius'])

    '''
    def _set_s2n_pdf(self):
        s2nspec = self['obj_model']['s2n']

        if s2nspec['type']=='uniform':
            s2n_r = s2nspec['range']
            self.s2n_pdf=ngmix.priors.FlatPrior(s2n_r[0], s2n_r[1])
        elif s2nspec['type']=='lognormal':
            self.s2n_pdf=ngmix.priors.LogNormal(s2nspec['mean'],s2nspec['sigma'])
        else:
            raise ValueError("bad s2n pdf type: '%s'" % s2nspec['type'])
    '''

    def _set_flux_pdf(self):
        fluxspec = self['obj_model']['flux']

        self.flux_is_in_log = fluxspec.get('is_in_log',False)
        if self.flux_is_in_log:
            print("Flux pdf is log")

        if fluxspec['type']=='uniform':
            flux_r = fluxspec['range']
            self.flux_pdf=ngmix.priors.FlatPrior(flux_r[0], flux_r[1])
        elif fluxspec['type']=='lognormal':
            self.flux_pdf=ngmix.priors.LogNormal(fluxspec['mean'],fluxspec['sigma'])
        elif fluxspec['type']=='gmixnd':
            self.flux_pdf=load_gmixnd(fluxspec)
        else:
            raise ValueError("bad flux pdf type: '%s'" % fluxspec['type'])

class SimBD(SimGS):
    """
    specific sim to deal with complications of a bulge+disk model
    """
    def _get_gal_obj(self):
        """
        get the galsim object for the galaxy model
        """

        pars=self._get_galaxy_pars()

        flux = pars['flux']

        r50 = pars['r50']

        g1,g2=pars['g']
        shear, shindex = self.shear_pdf.get_shear()
        pars['shear'] = shear
        pars['shear_index'] = shindex


        cenoff = pars['cenoff']

        fracdev = pars['fracdev']
        disk_flux = flux*(1.0 - fracdev)
        bulge_flux = flux*fracdev

        disk = galsim.Exponential(flux=disk_flux, half_light_radius=r50)
        bulge = galsim.DeVaucouleurs(flux=bulge_flux, half_light_radius=r50)

        # both disk and bulge get same overall shape
        disk  = disk.shear(g1=g1, g2=g2)
        bulge = bulge.shear(g1=g1, g2=g2)

        # the bulge can be offset from the disk
        dev_offset = pars['dev_offset']
        if dev_offset is not None:
            bulge = bulge.shift(dx=dev_offset[0], dy=dev_offset[1])

        # combine them and shear that
        gal = galsim.Add([disk, bulge])

        gal = gal.shear(g1=shear.g1, g2=shear.g2)

        if cenoff is not None:
            gal = gal.shift(dx=cenoff[0], dy=cenoff[1])

        tup=(r50,fracdev,dev_offset,cenoff)
        print("    r50: %g fracdev: %g dev_offset: %s cenoff: %s" % tup)

        return gal, pars

    def _get_galaxy_pars(self):
        """
        all pars are the same except for the shift of the bulge
        """
        pars=super(SimBD,self)._get_galaxy_pars()

        pars['fracdev'] = self.fracdev_pdf.sample()

        # distribution is in units of r50
        r50 = pars['r50']

        if self.dev_offset_pdf is not None:
            dev_offset1,dev_offset2 = self.dev_offset_pdf.sample2d()
            dev_offset = (r50*dev_offset1, r50*dev_offset2)
            print("    dev offset:",dev_offset[0],dev_offset[1])
        else:
            dev_offset=None
        pars['dev_offset'] = dev_offset

        return pars

    def _set_pdfs(self):
        """
        add fracdev and bulge offset distributions
        """
        super(SimBD,self)._set_pdfs()

        self._set_dev_offset_pdf()
        self._set_fracdev_pdf()

    def _set_dev_offset_pdf(self):
        ds_spec=self['obj_model']['dev_shift']
        if ds_spec is not None:
            # radius in units of r50
            self.dev_offset_pdf = ngmix.priors.ZDisk2D(ds_spec['radius'])
        else:
            self.dev_offset_pdf = None

    def _set_fracdev_pdf(self):
        bdr = self['obj_model']['fracdev']['range']
        self.fracdev_pdf=ngmix.priors.FlatPrior(bdr[0], bdr[1])


def quick_fit_gauss(image, maxiter=4000, tol=1.0e-6, ntry=4):
    """
    use EM to fit a single gaussian
    """

    dims = numpy.array(image.shape)
    cenguess = (dims-1)/2.0

    j=ngmix.UnitJacobian(0.0, 0.0)

    obs = ngmix.Observation(image, jacobian=j)

    guess_T = 4.0

    for i in xrange(ntry):
        guess_pars = [cenguess[0] + 0.1*srandu(),
                      cenguess[1] + 0.1*srandu(),
                      0.0 + 0.02*srandu(),
                      0.0 + 0.02*srandu(),
                      guess_T*(1.0 + 0.05*srandu()),
                      1.0 + 0.05*srandu() ]

        guess=ngmix.gmix.GMixModel(guess_pars, "gauss")

        fitter=ngmix.em.fit_em(obs, guess, maxiter=maxiter, tol=tol)

        res=fitter.get_result()
        if res['flags']==0:
            break
    
    if res['flags'] != 0:
        raise TryAgainError("could not fit 1 gauss")

    #pprint.pprint(res)
    return fitter
