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

from .util import TryAgainError

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
            impair = si.get_image_pair()
            # process images

        for i in xrange(1000):
            impair = si.get_image()
            # process images
        """

        self.update(sim_conf)

        seed = ngmixsim.get_devrand_uint()
        numpy.random.seed(seed)

        #base_rng = galsim.BaseDeviate(seed)
        #self.gauss_noise = galsim.GaussianNoise(base_rng, sigma=self['skysig'])
        #self.base_rng = base_rng

        self._setup()

        pprint.pprint(self)

    def get_image_pair(self):
        """
        Get an image pair, with noise added

        We never use a ring test for this version of the sim,
        so this is just for backwards compatibility
        """

        im1_dict=self.get_image()
        im2_dict=self.get_image()

        imdict={'im1':im1_dict, 'im2':im2_dict}

        return imdict

    def get_image(self):
        """
        get a randomized galaxy image
        """

        s2n = self.s2n_pdf.sample()
        print("    s2n: %g" % s2n)

        gal_pars=self._get_galaxy_pars()

        gal, psf = self._get_galsim_objects(gal_pars)

        psf_obs = self._make_obs(psf, self['psf']['s2n'])
        gal_obs = self._make_obs(gal, s2n)

        gal_obs.set_psf(psf_obs)

        save_pars=[gal_pars['r50'],
                   gal_pars['fracdev'],
                   gal_obs.image_flux]
        return {'obs':gal_obs,
                's2n':s2n,
                'model':self['model'],
                'gal_info':gal_pars,
                'pars':save_pars}

    def _make_obs(self, gs_obj, s2n):
        """
        get an ngmix Observation
        """

        dims = self['stamp_size']
        gsimage = galsim.ImageD(dims[0], dims[1])

        gs_obj.drawImage(gsimage, scale=self['pixel_scale'])

        im0 = gsimage.array
        image_nonoise, image, flux = self._scale_and_add_noise(im0, s2n)

        jacob = self._get_jacobian(image_nonoise)

        weight = numpy.zeros( image.shape ) + self['ivar']

        obs = ngmix.Observation(image, weight=weight, jacobian=jacob)

        # monkey patching
        obs.image_nonoise = image_nonoise
        obs.image_flux = flux

        if False:
            self._compare_images(image_nonoise,image,label1='im',label2='noisy')

        return obs

    def _compare_images(self, im1, im2, **keys):
        import images
        keys['width']=1000
        keys['height']=1000
        images.compare_images(im1, im2, **keys)


    def _scale_and_add_noise(self, im0, s2n):
        """
        find the flux that gives the requested s2n.
        """

        flux = self._get_flux_from_s2n(im0, s2n)

        imsum = im0.sum()

        factor = flux/imsum

        scaled_image = im0 * factor

        nim = numpy.random.normal(loc=0.0,
                                  scale=self['skysig'],
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

        scale=self['pixel_scale']
        return ngmix.Jacobian(row, col,
                              scale,
                              0.0,
                              0.0,
                              scale)

    def _get_galsim_objects(self, gal_pars):
        """
        get the galaxy and psf galsim objects
        """

        psf  = self._get_psf_obj(gal_pars['cenoff'])
        gal0 = self._get_gal_obj(gal_pars)

        gal = galsim.Convolve([psf, gal0])

        return gal, psf

    def _get_gal_obj(self, pars):
        """
        get the galsim object for the galaxy model
        """


        r50 = pars['r50']
        fracdev = pars['fracdev']

        g1,g2=pars['g']
        s1,s2=self['shear']

        disk_flux = 1.0 - fracdev
        bulge_flux = fracdev

        disk = galsim.Exponential(flux=disk_flux, half_light_radius=r50)
        bulge = galsim.DeVaucouleurs(flux=bulge_flux, half_light_radius=r50)

        # both disk and bulge get same overall shape
        disk  = disk.shear(g1=g1, g2=g2)
        bulge = bulge.shear(g1=g1, g2=g2)

        # there can be a shift to the centroid
        cenoff = pars['cenoff']
        if cenoff is not None:
            disk = disk.shift(dx=cenoff[0], dy=cenoff[1])
            bulge = bulge.shift(dx=cenoff[0], dy=cenoff[1])

        # the bulge can be offset from the disk
        dev_offset = pars['dev_offset']
        if dev_offset is not None:
            bulge = bulge.shift(dx=dev_offset[0], dy=dev_offset[1])

        # combine them and shear that
        gal = galsim.Add([disk, bulge])
        gal = gal.shear(g1=s1, g2=s2)

        return gal

    def _get_psf_obj(self, cenoff):
        """
        get the galsim object for the psf
        """

        pspec = self['psf']

        model=pspec['model']
        if model=='moffat':
            psf = galsim.Moffat(beta=pspec['beta'],
                                half_light_radius=pspec['r50'])
        elif model=='gauss':
            psf = galsim.Gaussian(half_light_radius=pspec['r50'])
        else:
            raise ValueError("bad psf model: '%s'" % model)

        psf = psf.shear(g1=pspec['shape'][0],
                        g2=pspec['shape'][1])

        if self.cen_pdf is not None:
            psf = psf.shift(dx=cenoff[0], dy=cenoff[1])
        return psf

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

        fracdev = self.fracdev_pdf.sample()

        # distribution is in units of r50
        if self.dev_offset_pdf is not None:
            dev_offset1,dev_offset2 = self.dev_offset_pdf.sample2d()
            dev_offset = (r50*dev_offset1, r50*dev_offset2)
            print("    dev offset:",dev_offset[0],dev_offset[1])
        else:
            dev_offset=None

        pars = {'model':self['model'],
                'g':(g1,g2),
                'r50':r50,
                'fracdev':fracdev,
                'cenoff':cenoff,
                'dev_offset':dev_offset}

        #pprint.pprint(pars)

        return pars

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
        self['pixel_scale'] = self.get('pixel_scale',1.0)
        self['ivar'] = 1.0/self['skysig']**2
        self['model'] = self['obj_model']['model']
        self._set_pdfs()

    def _set_pdfs(self):
        """
        Set all the priors
        """

        omodel=self['obj_model']

        s2n_r = omodel['s2n']['range']
        self.s2n_pdf=ngmix.priors.FlatPrior(s2n_r[0], s2n_r[1])

        g_spec=omodel['g']
        self.g_pdf=ngmix.priors.GPriorBA(g_spec['sigma'])

        r50_r = omodel['r50']['range']
        self.r50_pdf=ngmix.priors.FlatPrior(r50_r[0], r50_r[1])

        cr=omodel['cen_shift']
        if cr is None:
            self.cen_pdf=None
        else:
            self.cen_pdf=ngmix.priors.FlatPrior(-cr['radius'], cr['radius'])

        ds_spec=omodel['dev_shift']
        if ds_spec is not None:
            # radius in units of r50
            self.dev_offset_pdf = ngmix.priors.ZDisk2D(ds_spec['radius'])
        else:
            self.dev_offset_pdf = None

        bdr = omodel['fracdev']['range']
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
