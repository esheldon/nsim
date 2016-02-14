"""
Simulate images using galsim instead of ngmix
"""
from __future__ import print_function
import os
from pprint import pprint

import numpy

import ngmix
from . import sim as ngmixsim
from .sim import srandu

from .util import TryAgainError, load_gmixnd

from .pdfs import DiscreteSampler

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

        for k in self:
            if k != "shear":
                pprint(self[k])

        if 'seed' in self:
            print("    using seed:",self['seed'])
            numpy.random.seed(self['seed'])


    def get_image(self):
        """
        get a randomized galaxy image
        """

        gal, gal_pars, psf, psf_pars = self._get_galsim_objects()

        if 'psf_stamp_size' in self:
            nrows,ncols=self['psf_stamp_size']
        else:
            nrows,ncols=self['stamp_size']
        psf_obs = self._make_obs(psf, nrows, ncols, s2n=self['psf']['s2n'], isgal=False)

        nrows,ncols=self['stamp_size']
        gal_obs = self._make_obs(gal, nrows, ncols, isgal=True)

        s2n = self._get_expected_s2n(gal_obs.image_nonoise)
        print("    s2n expected:",s2n)

        gal_obs.set_psf(psf_obs)

        save_pars=[gal_pars['size'], gal_pars['flux']]

        if psf_pars['fwhm'] is None:
            psf_save_pars=psf_pars['r50']
        else:
            psf_save_pars=psf_pars['fwhm']

        return {'obs':gal_obs,
                's2n':s2n,
                'model':self['model'],
                'gal_info':gal_pars,
                'pars':save_pars,
                'psf_pars':psf_save_pars,
                'psf_obj': psf,
                'gal_obj': gal}


    def _make_obs(self, gs_obj, nrows, ncols, s2n=None, isgal=True):
        """
        get an ngmix Observation

        for psfs, send s2n=
        """

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

        bmask=None
        if isgal and self['bad_pixels'] is not None:
            bmask=self._set_bad_pixels(image, weight)
        elif isgal and self['masks'] is not None:
            bmask=self._set_mask(image, weight)

        obs = ngmix.Observation(
            image,
            weight=weight,
            bmask=bmask,
            jacobian=jacob
        )

        # monkey patching
        obs.image_nonoise = image_nonoise

        if False and s2n is None:
            self._compare_images(image_nonoise,image,label1='im',label2='noisy')

        return obs

    def _set_bad_pixels(self, image, weight):
        """
        currently one per stamp, random location
        """
        bmask=numpy.zeros(image.shape, dtype='i2')

        # ravel returns a view
        imravel = image.ravel()
        wtravel = weight.ravel()
        bmravel = bmask.ravel()

        ibad = numpy.random.randint(0, imravel.size)

        if self['bad_pixels']['replace_with']=='noise':
            imravel[ibad] = numpy.random.normal(
                loc=0.0,
                scale=self['noise'],
                size=1
            )
        else:
            imravel[ibad] = -9999.0

        wtravel[ibad] = 0.0
        bmravel[ibad] = 1

        return bmask

    def _set_mask(self, image, weight):
        """
        currently one per stamp, random location
        """
        mask = self._get_random_mask()

        bmask = numpy.zeros(mask.shape, dtype='i2')

        # ravel returns a view
        #imorig=image.copy()
        imravel = image.ravel()
        wtravel = weight.ravel()

        mask_ravel = mask.ravel()
        bmravel = bmask.ravel()

        mess="mask size does not match image"
        assert imravel.size == mask_ravel.size,mess

        ibad, = numpy.where(mask_ravel == 0)

        rep=self['masks']['replace_with']
        if rep =='noise':
            imravel[ibad] = numpy.random.normal(
                loc=0.0,
                scale=self['noise'],
                size=ibad.size
            )
        else:
            imravel[ibad] = rep

        wtravel[ibad] = 0.0
        bmravel[ibad] = 1

        if False:
            import images
            images.view_mosaic([mask,image,image-imorig],
                               titles=['mask','mod image','mod-orig'])
            if raw_input('hit a key: ') == 'q':
                stop

        return bmask


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

        gal_pars=self._get_galaxy_pars()

        psf, psf_pars  = self._get_psf_obj(gal_pars['cenoff'])

        if gal_pars['model']=='star':
            gal = psf.withFlux(gal_pars['flux'])
        else:
            gal0 = self._get_gal_obj(gal_pars)
            gal = galsim.Convolve([psf, gal0])

        return gal, gal_pars, psf, psf_pars

    def _get_gal_obj(self, pars):
        """
        get the galsim object for the galaxy model
        """

        flux = pars['flux']

        r50 = pars['r50']

        g1,g2=pars['g']

        cenoff = pars['cenoff']

        if pars['model']=='gauss':
            gal = galsim.Gaussian(flux=flux, half_light_radius=r50)
        elif pars['model']=='exp':
            gal = galsim.Exponential(flux=flux, half_light_radius=r50)
        elif pars['model']=='dev':
            gal = galsim.DeVaucouleurs(flux=flux, half_light_radius=r50)
        else:
            raise ValueError("bad galaxy model: '%s'" % pars['model'])

        # first give it an intrinsic shape
        gal = gal.shear(g1=g1, g2=g2)

        # now shear it
        if 'shear' in pars:
            shear=pars['shear']
            gal = gal.shear(g1=shear.g1, g2=shear.g2)

        # in the demos, the shift was always applied after the shear, not sure
        # if it matters
        if cenoff is not None:
            gal = gal.shift(dx=cenoff[0], dy=cenoff[1])

        tup=(r50,cenoff)
        print("    r50: %g cenoff: %s" % tup)

        return gal


    def _get_psf_obj(self, cenoff):
        """
        get the galsim object for the psf
        """

        pspec = self['psf']

        model=pspec['model']

        r50, fwhm = self._get_psf_size()

        if model=='moffat':
            psf = galsim.Moffat(beta=pspec['beta'],
                                half_light_radius=r50,
                                fwhm=fwhm)
        elif model=='gauss':
            psf = galsim.Gaussian(half_light_radius=r50,
                                  fwhm=fwhm)
        else:
            raise ValueError("bad psf model: '%s'" % model)

        psf_g1, psf_g2 = self._get_psf_shape()
        psf = psf.shear(g1=psf_g1, g2=psf_g2)

        if cenoff is not None:
            psf = psf.shift(dx=cenoff[0], dy=cenoff[1])

        return psf, {'fwhm':fwhm, 'r50':r50}

    def _get_psf_size(self):
        r50=None
        fwhm=None
        if self.psf_r50_pdf is not None:
            r50 = self.psf_r50_pdf.sample()
            print("    psf r50: %g" % r50)
        elif self.psf_fwhm_pdf is not None:
            fwhm = self.psf_fwhm_pdf.sample()
            print("    psf fwhm: %g" % fwhm)
        else:
            r50 = self['psf']['r50']
        return r50, fwhm

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

        if self.g_pdf is not None:
            g1,g2 = self.g_pdf.sample2d()
        else:
            g1,g2=None,None

        flux = self.flux_pdf.sample()
        if self.flux_is_in_log:
            flux = numpy.exp(flux)

        # this is the round r50
        if self.r50_pdf is not None:
            r50 = self.r50_pdf.sample()
        else:
            r50=None


        pars = {'model':self['model'],
                'g':(g1,g2),
                'flux':flux,
                'r50':r50,
                'size':r50,
                'cenoff':cenoff}

        if self.shear_pdf is not None:
            shear,shindex = self.shear_pdf.get_shear()
            pars['shear'] = shear
            pars['shear_index'] = shindex

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

        self['bad_pixels'] = self.get('bad_pixels',None)
        if self['bad_pixels'] is not None:
                assert self['bad_pixels']['replace_with']=='noise'

        self['masks'] = self.get('masks',None)
        if self['masks'] is not None:
            self._load_masks()

        self._set_pdfs()

    def _load_masks(self):
        """
        load masks from a multi-extension fits file
        """
        import fitsio
        mask_file=self['masks']['mask_file']

        print("Loading masks from:",mask_file)
        with fitsio.FITS(mask_file) as fits:
            mask_list = [hdu.read() for hdu in fits]

        self.mask_list=mask_list

    def _get_random_mask(self):
        """
        get a random mask, rotated by a random multiple
        of 90 degrees
        """
        i = numpy.random.randint(0, len(self.mask_list))
        #print("    loading mask:",i)
        mask=self.mask_list[i]

        ri=numpy.random.randint(0,4)
        mask = numpy.rot90(mask, k=ri)

        return mask

    def _set_pdfs(self):
        """
        Set all the priors
        """

        self._set_psf_pdf()
        self._set_flux_pdf()
        self._set_g_pdf()
        self._set_size_pdf()
        self._set_cen_pdf()
        self._set_shear_pdf()

    def _set_shear_pdf(self):
        from .shearpdf import ConstShearGenerator

        if 'shear' in self:
            shconf = self['shear']
            if shconf['type'] == 'const':
                pdf = ConstShearGenerator(shconf['shears'])
            else:
                raise ValueError("only shear 'const' for now")

            self.shear_pdf=pdf
        else:
            self.shear_pdf=None

    def _set_psf_pdf(self):
        pspec = self['psf']

        self.psf_r50_pdf=None
        self.psf_fwhm_pdf=None

        if 'fwhm' in pspec:
            assert pspec['fwhm']['type']=='discrete-pdf'
            fname=os.path.expandvars( pspec['fwhm']['file'] )
            print("Reading fwhm values from file:",fname)
            vals=numpy.fromfile(fname, sep='\n')
            self.psf_fwhm_pdf=DiscreteSampler(vals)
        else:
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
        if 'g' in self['obj_model']:
            g_spec=self['obj_model']['g']
            self.g_pdf=ngmix.priors.GPriorBA(g_spec['sigma'])
        else:
            self.g_pdf=None

    def _set_size_pdf(self):
        if 'r50' in self['obj_model']:
            r50spec = self['obj_model']['r50']

            if r50spec['type']=='uniform':
                r50_r = r50spec['range']
                self.r50_pdf=ngmix.priors.FlatPrior(r50_r[0], r50_r[1])
            elif r50spec['type']=='lognormal':
                self.r50_pdf=ngmix.priors.LogNormal(r50spec['mean'],
                                                    r50spec['sigma'])
            else:
                raise ValueError("bad r50 pdf type: '%s'" % r50spec['type'])
        else:
            self.r50_pdf=None

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

class SimGMix(SimGS):
    """
    the galaxy model is described by a gaussian mixture with size given
    by T instead of r50 etc.
    """

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

        if self.g_pdf is not None:
            g1,g2 = self.g_pdf.sample2d()
        else:
            g1,g2=None,None

        flux = self.flux_pdf.sample()
        if self.flux_is_in_log:
            flux = numpy.exp(flux)

        # this is the round T
        if self.T_pdf is not None:
            T = self.T_pdf.sample()
        else:
            T=None


        pars = {'model':self['model'],
                'g':(g1,g2),
                'flux':flux,
                'T':T,
                'size':T,
                'cenoff':cenoff}

        if self.shear_pdf is not None:
            shear,shindex = self.shear_pdf.get_shear()
            pars['shear'] = shear
            pars['shear_index'] = shindex

        return pars

    def _get_gal_obj(self, pars):
        """
        get the galsim object for the galaxy model
        """

        flux = pars['flux']

        T = pars['T']

        g1,g2=pars['g']

        cenoff = pars['cenoff']

        gm_pars=[0.0, 0.0, 0.0, 0.0, T, flux]
        gm = ngmix.GMixModel(gm_pars, self['obj_model']['gmix_model']) 

        gal = gm.make_galsim_object()

        # first give it an intrinsic shape
        gal = gal.shear(g1=g1, g2=g2)

        # now shear it
        if 'shear' in pars:
            shear=pars['shear']
            gal = gal.shear(g1=shear.g1, g2=shear.g2)

        # in the demos, the shift was always applied after the shear, not sure
        # if it matters
        if cenoff is not None:
            gal = gal.shift(dx=cenoff[0], dy=cenoff[1])

        tup=(T,cenoff)
        print("    T: %g cenoff: %s" % tup)

        return gal


    def _set_size_pdf(self):
        assert 'T' in self['obj_model']

        spec = self['obj_model']['T']

        if spec['type']=='lognormal':
            self.T_pdf=ngmix.priors.LogNormal(spec['mean'],
                                              spec['sigma'])
        elif spec['type']=='gmixnd':
            self.T_pdf=load_gmixnd(spec)
        else:
            raise ValueError("bad r50 pdf type: '%s'" % r50spec['type'])


class SimBD(SimGS):
    """
    specific sim to deal with complications of a bulge+disk model
    """
    def _get_gal_obj(self, pars):
        """
        get the galsim object for the galaxy model
        """

        flux = pars['flux']

        r50 = pars['r50']

        g1,g2=pars['g']

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

        if 'shear' in pars:
            shear=pars['shear']
            gal = gal.shear(g1=shear.g1, g2=shear.g2)

        if cenoff is not None:
            gal = gal.shift(dx=cenoff[0], dy=cenoff[1])

        tup=(r50,fracdev,dev_offset,cenoff)
        print("    r50: %g fracdev: %g dev_offset: %s cenoff: %s" % tup)

        return gal

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
        else:
            dev_offset=None
        pars['dev_offset'] = dev_offset

        if self.shear_pdf is not None:
            shear, shindex = self.shear_pdf.get_shear()
            pars['shear'] = shear
            pars['shear_index'] = shindex

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

class SimBDD(SimBD):
    """
    different ellipticities
    """
    def _get_gal_obj(self, pars):
        """
        get the galsim object for the galaxy model
        """

        flux = pars['flux']

        r50 = pars['r50']

        disk_g=pars['disk_g']
        bulge_g=pars['bulge_g']

        cenoff = pars['cenoff']

        fracdev = pars['fracdev']
        disk_flux = flux*(1.0 - fracdev)
        bulge_flux = flux*fracdev

        disk = galsim.Exponential(flux=disk_flux, half_light_radius=r50)
        bulge = galsim.DeVaucouleurs(flux=bulge_flux, half_light_radius=r50)

        # both disk and bulge get same overall shape
        disk  = disk.shear(g1=disk_g[0], g2=disk_g[1])
        bulge = bulge.shear(g1=bulge_g[0], g2=bulge_g[1])

        # the bulge can be offset from the disk
        dev_offset = pars['dev_offset']
        if dev_offset is not None:
            bulge = bulge.shift(dx=dev_offset[0], dy=dev_offset[1])

        # combine them and shear that
        gal = galsim.Add([disk, bulge])

        if 'shear' in pars:
            shear=pars['shear']
            gal = gal.shear(g1=shear.g1, g2=shear.g2)

        if cenoff is not None:
            gal = gal.shift(dx=cenoff[0], dy=cenoff[1])

        tup=(r50,fracdev,dev_offset,cenoff)
        print("    r50: %g fracdev: %g dev_offset: %s cenoff: %s" % tup)
        print("    disk g:",disk_g,"bulge g:",bulge_g)

        return gal

    def _get_galaxy_pars(self):
        """
        all pars are the same except for the shift of the bulge
        """

        # this will fill in only one of the shapes
        pars=super(SimBDD,self)._get_galaxy_pars()

        pars['disk_g'] = pars['g']

        if self.g_pdf is not None:
            angle =self.bulge_rot_pdf.sample()
            frac=self['obj_model']['bulge_gfrac']
            sh = ngmix.Shape(pars['disk_g'][0]*frac,
                             pars['disk_g'][1]*frac)
            shrot=sh.get_rotated(angle)
            pars['bulge_g'] = (shrot.g1, shrot.g2)
        else:
            pars['bulge_g'] = None

        return pars

    def _set_pdfs(self):
        """
        add fracdev and bulge offset distributions
        """
        super(SimBDD,self)._set_pdfs()

        self._set_bulge_rot_pdf()

    def _set_bulge_rot_pdf(self):
        sigma=numpy.deg2rad( self['obj_model']['bulge_rot_sigma_degrees'] )
        self.bulge_rot_pdf = ngmix.priors.Normal(0.0, sigma)


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

    return fitter
