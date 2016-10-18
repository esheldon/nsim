"""
Simulate images using galsim instead of ngmix
"""
from __future__ import print_function
import os, sys
from pprint import pprint

import numpy
import fitsio

import ngmix
import esutil as eu
from esutil.numpy_util import between

from . import sim as ngmixsim
from ngmix.priors import srandu

from .util import TryAgainError, load_gmixnd

from . import pdfs
from .pdfs import DiscreteSampler, PowerLaw, DiscreteHLRFluxSampler
import galsim
#try:
#    import galsim
#    have_galsim=True
#except ImportError:
#    have_galsim=False

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

        seed=self.get('seed',None)
        print("    using seed:",self['seed'])

        # seeding both the global and the local rng.  With the
        # local, we produce the same sim independent of the fitting
        # code which may use the global.
        numpy.random.seed(seed)
        self.rng=numpy.random.RandomState(seed=numpy.random.randint(0,2**30))
        self.galsim_rng = galsim.BaseDeviate(numpy.random.randint(0,2**30))

        self._setup()

        for k in self:
            if k != "shear":
                pprint(self[k])

        self.counter=0

    def get_image(self):
        """
        get a randomized galaxy image
        """

        wcs=self.get_galsim_wcs()

        gal, gal_pars, psf, psf_pars = self._get_galsim_objects()

        if 'stamp_size' in self:
            nrows,ncols=self['stamp_size']
        else:
            nrows,ncols=None,None

        psf_nrows,psf_ncols=self['psf_stamp_size']

        cen_shift=self._get_cen_shift()

        psf_obs, pflux = self._make_obs(psf, psf_nrows, psf_ncols, wcs, cen_shift,
                                        s2n=self['psf']['s2n'], isgal=False)

        # this will be none for sims where we generate a flux pdf
        s2n=gal_pars['s2n']
        gal_obs, gflux = self._make_obs(gal, nrows, ncols, wcs, cen_shift,
                                        s2n=s2n, isgal=True)

        if gflux is not None:
            gal_pars['flux'] = gflux

        s2n = self._get_expected_s2n(gal_obs.image_nonoise)
        print("    s2n expected:",s2n)

        gal_obs.set_psf(psf_obs)

        save_pars=[
            gal_pars['size'],
            gal_pars['flux'],
            #gal_pars['g'][0],
            #gal_pars['g'][1],
        ]

        if psf_pars['fwhm'] is None:
            psf_save_pars=psf_pars['r50']
        else:
            psf_save_pars=psf_pars['fwhm']

        self.counter += 1
        return {'obs':gal_obs,
                's2n':s2n,
                'model':self['model'],
                'gal_info':gal_pars,
                'pars':save_pars,
                'psf_pars':psf_save_pars,
                'psf_obj': psf,
                'gal_obj': gal}


    def _get_cen_shift(self):
        if self.cen_pdf is not None:
            coff1 = self.cen_pdf.sample()
            coff2 = self.cen_pdf.sample()

            cen_shift=(coff1,coff2)
        else:
            cen_shift=None

        if cen_shift is not None:
            print("    cen_shift: %g,%g" % cen_shift)

        return cen_shift

    def _make_obs(self, gs_obj, nrows, ncols, wcs, cen_shift, s2n=None, isgal=True):
        """
        get an ngmix Observation

        for psfs, send s2n=
        """

        #gsimage = galsim.ImageD(ncols, nrows, wcs=self.galsim_wcs)
        #gs_obj.drawImage(image=gsimage)
        gsimage = gs_obj.drawImage(nx=ncols,
                                   ny=nrows,
                                   wcs=wcs,
                                   dtype=numpy.float64,
                                   offset=cen_shift)
        #print("    dims: %d,%d" % tuple(gsimage.array.shape))
        im0 = gsimage.array
        if s2n is not None:
            image_nonoise, image, flux = self._scale_and_add_noise(im0, s2n)
        else:
            flux=None
            image_nonoise = im0.copy()
            image = self._add_noise(image_nonoise)

        jacob = self._get_jacobian(image_nonoise,gsimage.wcs)

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

        return obs, flux

    def _set_bad_pixels(self, image, weight):

        bp=self['bad_pixels']
        bcrate = bp['bad_column_rate']
        bprate = bp['bad_pixel_rate']

        bsum = bcrate + bprate

        r = self.rng.uniform()
        if r < bcrate:
            print("        doing bad column")
            return self._set_bad_column(image, weight)
        elif r < bsum:
            print("        doing bad pixel")
            return self._set_bad_pixel(image, weight)

    def _set_bad_pixel(self, image, weight):
        """
        currently one per stamp, random location
        """
        bmask=numpy.zeros(image.shape, dtype='i2')

        # ravel returns a view
        imravel = image.ravel()
        wtravel = weight.ravel()
        bmravel = bmask.ravel()

        ibad = self.rng.randint(0, imravel.size)

        badval=self['bad_pixels']['replace_with']
        imravel[ibad] = badval

        wtravel[ibad] = 0.0
        bmravel[ibad] = 1

        return bmask

    def _set_bad_column(self, image, weight):
        """
        one per stamp, random location but not
        within 5 pixels of the center
        """

        cencol = int( (image.shape[1]-1.0)/2.0 )
        lowlim = cencol-4.0
        highlim = cencol+4.0

        bmask=numpy.zeros(image.shape, dtype='i2')

        while True:
            badcol=self.rng.randint(0,image.shape[1])

            if badcol < lowlim or badcol > highlim:
                break

        badval=self['bad_pixels']['replace_with']
        image[:,badcol] = badval
        weight[:,badcol] = 0.0
        bmask[:,badcol] = 1

        return bmask


    def _set_mask(self, image, weight):
        """
        currently one per stamp, random location
        """
        mask = self._get_random_mask()

        # make a local copy
        bmask = numpy.zeros(mask.shape, dtype='i2')
        bmask[:,:] = mask[:,:]

        # ravel returns a view
        imravel = image.ravel()
        wtravel = weight.ravel()
        bmravel = bmask.ravel()

        mess="mask size does not match image"
        assert imravel.size == bmravel.size,mess

        ibad, = numpy.where(bmravel != 0)

        rep=self['masks']['replace_with']
        if rep =='noise':
            imravel[ibad] = self.rng.normal(
                loc=0.0,
                scale=self['noise'],
                size=ibad.size
            )
        else:
            imravel[ibad] = rep

        wtravel[ibad] = 0.0

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
        images.compare_images(im1, im2, 
                              #file='/u/ki/esheldon/public_html/tmp/plots/tmp.png',
                              **keys)
        key=raw_input('hit a key: ')
        if key=='q':
            stop


    def _add_noise(self, im0):
        #nim = numpy.random.normal(loc=0.0,
        nim = self.rng.normal(loc=0.0,
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

        #nim = numpy.random.normal(loc=0.0,
        nim = self.rng.normal(loc=0.0,
                              scale=self['noise'],
                              size=im0.shape)
        noisy_image = scaled_image + nim

        return scaled_image, noisy_image, flux

    def get_galsim_wcs(self):
        if 'wcs' in self:
            ws=self['wcs']
            dudx=ws['dudx']
            dudy=ws['dudy']
            dvdx=ws['dvdx']
            dvdy=ws['dvdy']

            if 'dudx_std' in ws:
                rng=self.rng
                dudx = dudx + rng.normal(scale=ws['dudx_std'])
                dudy = dudy + rng.normal(scale=ws['dudy_std'])
                dvdx = dvdx + rng.normal(scale=ws['dvdx_std'])
                dvdy = dvdy + rng.normal(scale=ws['dvdy_std'])
        else:
            dudx=1.0
            dudy=0.0
            dvdx=0.0
            dvdy=1.0

        wcs=galsim.JacobianWCS(
            dudx,
            dudy,
            dvdx,
            dvdy,
        )


        # our sims are always in pixels for now
        """        
        det = numpy.abs( dudy*dvdx-dudx*dvdy )
        scale = numpy.sqrt(det)

        wcs=galsim.JacobianWCS(dudx/scale,
                               dudy/scale,
                               dvdx/scale,
                               dvdy/scale)
        """
        return wcs

    def _get_jacobian(self, image, wcs):
        """
        find the best center and set the jacobian center there
        """
        fitter = quick_fit_gauss(image, self.rng)
        row,col = fitter.get_gmix().get_cen()

        return ngmix.Jacobian(wcs=wcs,
                              row=row,
                              col=col)

    def _get_galsim_objects(self):
        """
        get the galaxy and psf galsim objects
        """

        gal_pars=self._get_galaxy_pars()

        psf, psf_pars  = self._get_psf_obj()

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
        if flux is None:
            # we will scale for s/n later
            flux=1.0

        r50 = pars['r50']

        g1,g2=pars['g']

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

        print("    r50: %g" % r50)

        return gal


    def _get_psf_obj(self):
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
        elif 'r50' in self['psf']:
            r50 = self['psf']['r50']
        elif 'fwhm' in self['psf']:
            fwhm = self['psf']['fwhm']
        return r50, fwhm

    def _get_psf_shape(self):
        pspec = self['psf']

        if self.psf_ellip_pdf is not None:
            psf_g1, psf_g2 = self.psf_ellip_pdf.sample()
            print("    psf (pdf) shape: %g %g" % (psf_g1, psf_g2))

        elif 'randomized_orientation' in pspec:
            ro=pspec['randomized_orientation']
            if ro["dist"]=="uniform":
                #angle = numpy.random.random()*2*numpy.pi
                angle = self.rng.uniform()*2*numpy.pi
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

        if self['do_ring'] and (self.counter % 2) != 0:
            print("doing ring")
            gold=self.old_pars['g']
            shape=ngmix.Shape(gold[0],gold[1])
            shape=shape.get_rotated(numpy.pi/2.0)
            g1,g2=shape.g1,shape.g2

            pars=self.old_pars
            pars['g'] = (g1,g2)

            return pars

        else:
            if self.g_pdf is not None:
                g1,g2 = self.g_pdf.sample2d()
            else:
                g1,g2=None,None

        if self.flux_pdf is not None:
            flux = self.flux_pdf.sample()
            if self.flux_is_in_log:
                flux = numpy.exp(flux)
        else:
            flux=None

        if self.s2n_pdf is not None:
            s2n = self.s2n_pdf.sample()
        else:
            s2n=None


        # this is the round r50
        if self.r50_pdf is not None:
            r50 = self.r50_pdf.sample()
        else:
            r50=None


        pars = {'model':self['model'],
                'g':(g1,g2),
                'flux':flux,
                's2n':s2n,
                'r50':r50,
                'size':r50}

        if self.shear_pdf is not None:
            shear,shindex = self.shear_pdf.get_shear(self.rng)
            pars['shear'] = shear
            pars['shear_index'] = shindex

        self.old_pars=pars
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

        self['do_ring'] = self.get('do_ring',False)

        self['ivar'] = 1.0/self['noise']**2
        self['model'] = self['obj_model']['model']

        self['bad_pixels'] = self.get('bad_pixels',None)

        self['masks'] = self.get('masks',None)
        if self['masks'] is not None:
            self._load_masks()

        self._set_pdfs()

    def _load_masks(self):
        """
        load masks from a multi-extension fits file

        if add_rotated is set, a 90 degree rotated version
        is added to cancel symmetries in the mask, such as
        bad columns
        """
        mask_file=self['masks']['mask_file']
        add_rotated=self['masks']['add_rotated']

        print("Loading masks from:",mask_file)

        mask_list=[]
        with fitsio.FITS(mask_file) as fits:

            for hdu in fits:
                mask = hdu.read()
                if add_rotated:
                    rm = numpy.rot90(mask)
                    mask = mask + rm

                mask_list.append( mask )

        self.mask_list=mask_list

    def _get_random_mask(self):
        """
        get a random mask, rotated by a random multiple
        of 90 degrees
        """
        #i = numpy.random.randint(0, len(self.mask_list))
        i = self.rng.randint(0, len(self.mask_list))
        #print("    loading mask:",i)
        mask=self.mask_list[i]

        #ri=self.rng.randint(0,4)
        #mask = numpy.rot90(mask, k=ri)

        return mask

    def _set_pdfs(self):
        """
        Set all the priors
        """

        self._set_psf_pdf()

        # should generally just either flux or s2n
        self._set_flux_pdf()
        self._set_s2n_pdf()

        self._set_g_pdf()
        self._set_size_pdf()
        self._set_cen_pdf()
        self._set_shear_pdf()

    def _set_shear_pdf(self):
        from . import shearpdf

        # self.rng only used for a Selector
        self.shear_pdf = shearpdf.get_shear_pdf(self)

    def _set_psf_pdf(self):
        pspec = self['psf']

        self.psf_r50_pdf=None
        self.psf_fwhm_pdf=None

        if 'fwhm' in pspec:
            if isinstance(pspec['fwhm'], dict):
                type=pspec['fwhm']['type']
                if type=='discrete-pdf':
                    fname=os.path.expandvars( pspec['fwhm']['file'] )
                    print("Reading fwhm values from file:",fname)
                    vals=numpy.fromfile(fname, sep='\n')
                    psf_fwhm_pdf=DiscreteSampler(vals,
                                                 rng=self.rng)
                elif type=='lognormal':
                    psf_fwhm_pdf=ngmix.priors.LogNormal(
                        pspec['fwhm']['mean'],
                        pspec['fwhm']['sigma'],
                        rng=self.rng,
                    )
                else:
                    raise ValueError("bad fwhm type: '%s'" % type)

                if 'range' in pspec['fwhm']:
                    bounds=pspec['fwhm']['range']
                    print("   bounding psf fwhm to:",bounds)
                    psf_fwhm_pdf = ngmix.priors.Bounded1D(psf_fwhm_pdf,bounds)

                self.psf_fwhm_pdf  = psf_fwhm_pdf 
        else:
            if isinstance(pspec['r50'], dict):
                r50pdf = pspec['r50']
                assert r50pdf['type']=="lognormal","r50 pdf log normal for now"

                self.psf_r50_pdf = ngmix.priors.LogNormal(
                    r50pdf['mean'],
                    r50pdf['sigma'],
                    rng=self.rng,
                )

        if isinstance(pspec['shape'],dict):
            ppdf=pspec['shape']
            assert ppdf['type']=="normal2d"
            self.psf_ellip_pdf=ngmix.priors.SimpleGauss2D(
                ppdf['cen'][0],
                ppdf['cen'][1],
                ppdf['sigma'][0],
                ppdf['sigma'][1],
                rng=self.rng,
            )
        else:
            self.psf_ellip_pdf=None

    def _set_g_pdf(self):
        if 'g' in self['obj_model']:
            g_spec=self['obj_model']['g']
            self.g_pdf=ngmix.priors.GPriorBA(
                g_spec['sigma'],
                rng=self.rng,
            )
        else:
            self.g_pdf=None

    def _set_size_pdf(self):
        if 'r50' in self['obj_model']:


            r50spec = self['obj_model']['r50']
            if not isinstance(r50spec,dict):
                self.r50_pdf=DiscreteSampler([r50spec], rng=self.rng)
            else:

                if r50spec['type']=='uniform':
                    r50_r = r50spec['range']
                    self.r50_pdf=ngmix.priors.FlatPrior(
                        r50_r[0],
                        r50_r[1],
                        rng=self.rng,
                    )
                elif r50spec['type']=='lognormal':
                    self.r50_pdf=ngmix.priors.LogNormal(
                        r50spec['mean'],
                        r50spec['sigma'],
                        rng=self.rng,
                    )
                elif r50spec['type']=='discrete-pdf':
                    fname=os.path.expandvars( r50spec['file'] )
                    print("Reading r50 values from file:",fname)
                    vals=fitsio.read(fname)
                    self.r50_pdf=DiscreteSampler(vals, rng=self.rng)

                else:
                    raise ValueError("bad r50 pdf type: '%s'" % r50spec['type'])
        else:
            self.r50_pdf=None

    def _set_cen_pdf(self):

        cr=self['cen_shift']
        if cr is None:
            self.cen_pdf=None
        else:
            type=cr.get('type','uniform')
            if type=='uniform':
                self.cen_pdf=ngmix.priors.FlatPrior(
                    -cr['radius'], cr['radius'],
                    rng=self.rng,
                )
            else:
                raise ValueError("cen shift type should be 'uniform'")


    def _set_flux_pdf(self):

        if 'flux' not in self['obj_model']:
            self.flux_pdf=None
            return

        fluxspec = self['obj_model']['flux']
        if not isinstance(fluxspec,dict):
            self.flux_is_in_log = False
            self.flux_pdf=DiscreteSampler([fluxspec], rng=self.rng)
        else:


            self.flux_is_in_log = fluxspec.get('is_in_log',False)
            if self.flux_is_in_log:
                print("Flux pdf is log")

            if fluxspec['type']=='uniform':
                flux_r = fluxspec['range']
                self.flux_pdf=ngmix.priors.FlatPrior(
                    flux_r[0], flux_r[1],
                    rng=self.rng,
                )
            elif fluxspec['type']=='lognormal':
                self.flux_pdf=ngmix.priors.LogNormal(
                    fluxspec['mean'],
                    fluxspec['sigma'],
                    rng=self.rng,
                )
            elif fluxspec['type']=='gmixnd':
                self.flux_pdf=load_gmixnd(fluxspec,rng=self.rng)

            elif fluxspec['type']=='powerlaw':

                index=fluxspec['index']
                xmin=fluxspec['min']
                xmax=fluxspec['max']

                self.flux_pdf=PowerLaw(index, xmin, xmax)



            else:
                raise ValueError("bad flux pdf type: '%s'" % fluxspec['type'])

    def _set_s2n_pdf(self):

        if 's2n' not in self['obj_model']:
            self.s2n_pdf=None
            return

        s2nspec = self['obj_model']['s2n']
        if not isinstance(s2nspec,dict):
            self.s2n_pdf=DiscreteSampler([s2nspec], rng=self.rng)
        else:

            if s2nspec['type']=='uniform':
                s2n_r = s2nspec['range']
                self.s2n_pdf=ngmix.priors.FlatPrior(
                    s2n_r[0], s2n_r[1],
                    rng=self.rng,
                )
            elif s2nspec['type']=='lognormal':
                self.s2n_pdf=ngmix.priors.LogNormal(
                    s2nspec['mean'],
                    s2nspec['sigma'],
                    rng=self.rng,
                )
            elif s2nspec['type']=='gmixnd':
                self.s2n_pdf=load_gmixnd(s2nspec,rng=self.rng)

            elif s2nspec['type']=='powerlaw':

                index=s2nspec['index']
                xmin=s2nspec['min']
                xmax=s2nspec['max']

                self.s2n_pdf=PowerLaw(index, xmin, xmax)

            else:
                raise ValueError("bad s2n pdf type: '%s'" % s2nspec['type'])

    def _set_dev_offset_pdf(self):
        ds_spec=self['obj_model']['dev_shift']
        if ds_spec is not None:
            # radius in units of r50
            self.dev_offset_pdf = ngmix.priors.ZDisk2D(ds_spec['radius'],
                                                       rng=self.rng)
        else:
            self.dev_offset_pdf = None


class SimCosmos(SimGS):
    """
    using real galaxy fits from cosmos.
    """

    def _setup(self):

        self['ivar'] = 1.0/self['noise']**2
        self['model'] = self['obj_model']['model']
        assert self['obj_model']['type']=='parametric'

        self['bad_pixels'] = self.get('bad_pixels',None)

        self['masks'] = self.get('masks',None)
        if self['masks'] is not None:
            self._load_masks()

        self._set_pdfs()

        self._set_gal_cat()

    def _get_galaxy_pars(self):
        """
        all pars are the same except for the shift of the bulge
        """

        # compatibility with older code
        pars={
            'model':  self['model'],
            'flux':   -9999,
            'r50':    -9999,
            'size':   -9999,
            'g':      [-9999,-9999],
            's2n':    None,
        }


        if self.shear_pdf is not None:
            shear, shindex = self.shear_pdf.get_shear(self.rng)
            pars['shear'] = shear
            pars['shear_index'] = shindex


        return pars

    def _get_gal_obj(self, pars):
        """
        get the galsim object for the galaxy model
        """

        gal = self.cat.makeGalaxy(rng=self.galsim_rng)

        # there is a net orientation in these galaxies, we
        # need to rotate them

        theta = self.rng.uniform()*2.0*numpy.pi*galsim.radians
        gal = gal.rotate(theta)

        if 'shear' in pars:
            shear=pars['shear']
            gal = gal.shear(g1=shear.g1, g2=shear.g2)

        return gal



    def _set_pdfs(self):
        """
        Set all the priors
        """

        self._set_psf_pdf()

        # should generally just either flux or s2n
        self._set_cen_pdf()
        self._set_shear_pdf()

    def _set_gal_cat(self):
        self.cat=galsim.COSMOSCatalog(use_real=False)

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
                'size':T}

        if self.shear_pdf is not None:
            shear,shindex = self.shear_pdf.get_shear(self.rng)
            pars['shear'] = shear
            pars['shear_index'] = shindex

        return pars

    def _get_gal_obj(self, pars):
        """
        get the galsim object for the galaxy model
        """

        flux = pars['flux']
        if flux is None:
            # we will scale for s/n later
            flux=1.0

        T = pars['T']

        g1,g2=pars['g']

        gm_pars=[0.0, 0.0, 0.0, 0.0, T, flux]
        gm = ngmix.GMixModel(gm_pars, self['obj_model']['gmix_model']) 

        gal = gm.make_galsim_object()

        # first give it an intrinsic shape
        gal = gal.shear(g1=g1, g2=g2)

        # now shear it
        if 'shear' in pars:
            shear=pars['shear']
            gal = gal.shear(g1=shear.g1, g2=shear.g2)

        return gal


    def _set_size_pdf(self):
        assert 'T' in self['obj_model']

        spec = self['obj_model']['T']

        if spec['type']=='lognormal':
            self.T_pdf=ngmix.priors.LogNormal(
                spec['mean'],
                spec['sigma'],
                rng=self.rng,
            )
        elif spec['type']=='gmixnd':
            self.T_pdf=load_gmixnd(spec, rng=self.rng)
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
        if flux is None:
            # we will scale for s/n later
            flux=1.0

        r50 = pars['r50']

        g1,g2=pars['g']

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

        # always put the more prominent at the center
        if dev_offset is not None:
            bulge = bulge.shift(dx=dev_offset[0], dy=dev_offset[1])
        #    if fracdev > 0.5:
        #        disk = disk.shift(dx=dev_offset[0], dy=dev_offset[1])
        #    else:
        #        bulge = bulge.shift(dx=dev_offset[0], dy=dev_offset[1])

        # combine them and shear that
        gal = galsim.Add([disk, bulge])

        if 'shear' in pars:
            shear=pars['shear']
            gal = gal.shear(g1=shear.g1, g2=shear.g2)


        tup=(r50,fracdev,dev_offset)
        print("    r50: %g fracdev: %g dev_offset: %s" % tup)

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
            shear, shindex = self.shear_pdf.get_shear(self.rng)
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

    def _set_fracdev_pdf(self):
        bdr = self['obj_model']['fracdev']['range']
        self.fracdev_pdf=ngmix.priors.FlatPrior(bdr[0], bdr[1],
                                                rng=self.rng)



class SimBDJoint(SimBD):
    """
    joint flux-hlr distribution form cosmos sersic fits
    """
    def _set_pdfs(self):
        """
        Set all the priors
        """

        self._set_psf_pdf()

        # should generally just either flux or s2n
        self._set_joint_pdf()

        self._set_fracdev_pdf()
        self._set_dev_offset_pdf()
        self._set_g_pdf()
        self._set_cen_pdf()
        self._set_shear_pdf()

    def _get_galaxy_pars(self):
        """
        all pars are the same except for the shift of the bulge
        """

        s2n=None

        fr50 = self.joint_pdf.sample()
        r50=fr50['r50']
        flux=fr50['flux']

        fracdev = self.fracdev_pdf.sample()

        if self.g_pdf is not None:
            g1,g2 = self.g_pdf.sample2d()
        else:
            g1,g2=None,None

        if self.dev_offset_pdf is not None:
            dev_offset1,dev_offset2 = self.dev_offset_pdf.sample2d()
            dev_offset = (r50*dev_offset1, r50*dev_offset2)
        else:
            dev_offset=None

        pars = {
            'model':self['model'],
            'g':(g1,g2),
            'flux':flux,
            's2n':s2n,
            'r50':r50,
            'size':r50,
            'fracdev':fracdev,
            'dev_offset': dev_offset,
        }

        if self.shear_pdf is not None:
            shear,shindex = self.shear_pdf.get_shear(self.rng)
            pars['shear'] = shear
            pars['shear_index'] = shindex

        return pars


    def _get_cosmos_data(self):
        fname='real_galaxy_catalog_25.2_fits.fits'
        fname=os.path.join(
            sys.exec_prefix,
            'share',
            'galsim',
            'COSMOS_25.2_training_sample',
            fname,
        )

        fmin,fmax=self['obj_model']['flux_range']
        r50min,r50max=self['obj_model']['r50_range']

        print("reading cosmos file:",fname)
        data=fitsio.read(fname, lower=True)
        w,=numpy.where(
            (data['viable_sersic']==1) &
            #between(data['hlr'][:,0], 0.15, 3.0) &
            #between(data['flux'][:,0], 2.5, 100.0)
            between(data['hlr'][:,0], r50min, r50max) &
            between(data['flux'][:,0], fmin, fmax)
        )
        print("kept %d/%d" % (w.size, data.size))

        data=data[w]
        return data


    def _set_joint_pdf(self):
        """
        joint size-flux from the cosmos catalog
        """

        data=self._get_cosmos_data()
        self.joint_pdf=DiscreteHLRFluxSampler(data, rng=self.rng)


class SimBDJointDiffshape(SimBD):
    """
    joint flux-hlr distribution form cosmos sersic fits

    different ellipticities
    """
    def _set_pdfs(self):
        """
        Set all the priors
        """

        self._set_psf_pdf()

        # should generally just either flux or s2n
        self._set_joint_pdf()

        self._set_fracdev_pdf()
        self._set_g_pdf()
        self._set_cen_pdf()
        self._set_shear_pdf()

    def _set_g_pdf(self):
        gspec=self['obj_model']['gexp']
        self.gexp_pdf=ngmix.priors.GPriorBA(
            gspec['sigma'],
            rng=self.rng,
        )

        gspec=self['obj_model']['gdev']
        self.gdev_pdf=ngmix.priors.GPriorBA(
            gspec['sigma'],
            rng=self.rng,
        )

    def _get_gal_obj(self, pars):
        """
        get the galsim object for the galaxy model
        """

        flux = pars['flux']
        if flux is None:
            # we will scale for s/n later
            flux=1.0

        r50 = pars['r50']

        fracdev = pars['fracdev']
        disk_flux = flux*(1.0 - fracdev)
        bulge_flux = flux*fracdev

        disk = galsim.Exponential(flux=disk_flux, half_light_radius=r50)
        bulge = galsim.DeVaucouleurs(flux=bulge_flux, half_light_radius=r50)

        # both disk and bulge get same overall shape
        disk  = disk.shear(g1=pars['gexp'][0], g2=pars['gexp'][1])
        bulge = bulge.shear(g1=pars['gdev'][0], g2=pars['gdev'][1])

        olist=[disk,bulge]

        olist += self._make_knots(pars['knots'], r50)

        # combine them and shear that
        gal = galsim.Add(olist)

        if 'shear' in pars:
            shear=pars['shear']
            gal = gal.shear(g1=shear.g1, g2=shear.g2)


        tup=(r50,flux,fracdev)
        print("    r50: %g flux: %g fracdev: %g" % tup)

        return gal

    def _make_knots(self, knots, r50):
        olist=[]
        if knots is None:
            return olist

        for i in xrange(knots['num']):
            dx,dy = self.rng.normal(scale=r50, size=2)

            obj=galsim.Gaussian(sigma=1.0e-3, flux=knots['flux'])
            obj=obj.shift(dx=dx, dy=dy)

            olist.append(obj)

        return olist

    def _get_galaxy_pars(self):
        """
        all pars are the same except for the shift of the bulge
        """

        s2n=None

        r50, flux = self.joint_pdf.sample()
        fracdev = self.fracdev_pdf.sample()

        g1exp,g2exp = self.gexp_pdf.sample2d()
        g1dev,g2dev = self.gdev_pdf.sample2d()

        if 'knots' in self['obj_model']:
            kp=self['obj_model']['knots']
            num_knots=self.rng.poisson(lam=kp['nmean'])
            knot_flux = kp['flux_frac']*flux
            knots={'num':num_knots, 'flux':knot_flux}
            print("    num knots:",num_knots)
        else:
            knots=None

        pars = {
            'model':self['model'],
            'gexp':(g1exp,g2exp),
            'gdev':(g1dev,g2dev),
            'flux':flux,
            's2n':s2n,
            'r50':r50,
            'size':r50,
            'fracdev':fracdev,
            'knots':knots,
        }

        if self.shear_pdf is not None:
            shear,shindex = self.shear_pdf.get_shear(self.rng)
            pars['shear'] = shear
            pars['shear_index'] = shindex

        return pars

    def _set_joint_pdf(self):
        """
        joint size-flux from the cosmos catalog
        """

        self.joint_pdf = pdfs.CosmosR50Flux(
            self['obj_model']['r50_range'],
            self['obj_model']['flux_range'],
        )


    '''
    def _get_cosmos_data(self):
        fname='real_galaxy_catalog_25.2_fits.fits'
        fname=os.path.join(
            sys.exec_prefix,
            'share',
            'galsim',
            'COSMOS_25.2_training_sample',
            fname,
        )

        fmin,fmax=self['obj_model']['flux_range']
        r50min,r50max=self['obj_model']['r50_range']

        print("reading cosmos file:",fname)
        data=fitsio.read(fname, lower=True)
        w,=numpy.where(
            (data['viable_sersic']==1) &
            between(data['hlr'][:,0], r50min, r50max) &
            between(data['flux'][:,0], fmin, fmax)
        )
        #    between(data['hlr'][:,0], 0.15, 1.0) &
        #    between(data['flux'][:,0], 2.5, 100.0)
        print("kept %d/%d" % (w.size, data.size))

        data=data[w]
        return data

    def _set_joint_pdf(self):
        """
        joint size-flux from the cosmos catalog
        """

        data=self._get_cosmos_data()
        self.joint_pdf=DiscreteHLRFluxSampler(data, rng=self.rng)
    '''



class SimBDD(SimBD):
    """
    different ellipticities
    """
    def _get_gal_obj(self, pars):
        """
        get the galsim object for the galaxy model
        """

        flux = pars['flux']
        if flux is None:
            # we will scale for s/n later
            flux=1.0

        r50 = pars['r50']

        disk_g=pars['disk_g']
        bulge_g=pars['bulge_g']

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

        tup=(r50,fracdev,dev_offset)
        print("    r50: %g fracdev: %g dev_offset: %s" % tup)
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
        self.bulge_rot_pdf = ngmix.priors.Normal(0.0, sigma, rng=self.rng)


def quick_fit_gauss(image, rng, maxiter=4000, tol=1.0e-6, ntry=4):
    """
    use EM to fit a single gaussian
    """

    dims = numpy.array(image.shape)
    cenguess = (dims-1)/2.0

    j=ngmix.UnitJacobian(row=0.0, col=0.0)

    obs = ngmix.Observation(image, jacobian=j)

    guess_T = 4.0

    for i in xrange(ntry):
        guess_pars = [cenguess[0] + 0.1*srandu(rng=rng),
                      cenguess[1] + 0.1*srandu(rng=rng),
                      0.0 + 0.02*srandu(rng=rng),
                      0.0 + 0.02*srandu(rng=rng),
                      guess_T*(1.0 + 0.05*srandu(rng=rng)),
                      1.0 + 0.05*srandu(rng=rng) ]

        guess=ngmix.gmix.GMixModel(guess_pars, "gauss")

        fitter=ngmix.em.fit_em(obs, guess, maxiter=maxiter, tol=tol)

        res=fitter.get_result()
        if res['flags']==0:
            break

    if res['flags'] != 0:
        raise TryAgainError("could not fit 1 gauss")
    return fitter
