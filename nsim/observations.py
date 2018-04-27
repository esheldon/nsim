from __future__ import print_function
try:
    xrange
except:
    xrange=range
    raw_input=input

import os
import logging
import random
from copy import deepcopy
import numpy
import galsim
import ngmix
import fitsio
from . import pdfs

from .util import TryAgainError

logger = logging.getLogger(__name__)

BAD_PIXEL=2**0
BAD_COLUMN=2**1

def get_observation_maker(*args, **kw):
    conf = args[0]
    if 'nnbr' in conf and conf['nnbr'] > 0:
        #return NbrObservationMaker(*args)
        return NbrObservationMakerMulti(*args, **kw)
    else:
        return ObservationMaker(*args, **kw)
        

class ObservationMaker(dict):
    """
    images subset of the config

    noise: 1.0
    wcs:                    # optional
        dudx: 1.0
        dudy: 0.0
        dvdx: 0.0
        dvdy: 1.0

    psf:
        s2n: 10000

        stamp_size: [48,48] # optional

    object:
        nepoch: 10 # optional

        stamp_size: [48,48] # optional
        cen_shift:
            type: "uniform"
            radius: 0.5 # pixels
    """
    def __init__(self,
                 config,
                 psf_maker,
                 object_maker,
                 rng,
                 galsim_rng,
                 shear_pdf=None):

        self.update(config)

        self['use_canonical_center'] = self.get('use_canonical_center',False)

        self.psf_maker=psf_maker
        self.object_maker=object_maker
        self.shear_pdf=shear_pdf
        self.rng=rng
        self.galsim_rng=galsim_rng

        self._set_pdfs()
        self._load_bmasks()

    def __call__(self):

        objconf=self['object']
        nepoch = objconf.get('nepoch',1)

        cobjlist, psflist, wcslist, meta = self._get_convolved_object_info(
            objconf,
        )

        dims, psf_dims = self._get_dims(cobjlist, psflist, wcslist)

        obslist = ngmix.observation.ObsList()

        for i in xrange(nepoch):

            obs = self._get_obs(
                psflist[i],
                cobjlist[i],
                wcslist[i],
                dims,
                psf_dims,
            )

            obslist.append( obs )

        meta['s2n']  = get_expected_s2n(obslist)

        obslist.update_meta_data(meta)

        if 'defects' in self:
            self._add_defects(obslist)

        if 'coadd' in self:
            obslist = self._do_coadd(obslist)

        return obslist

    def _get_dims(self, cobjlist, psflist, wcslist):
        """
        get maximum size from all listed convolved objects
        """

        if self.nrows is None or self.psf_nrows is None:

            n=len(cobjlist)
            size=-1
            psf_size=-1
            for i in xrange(n):

                scale, shear, theta, flip = wcslist[i].getDecomposition()

                tsize=cobjlist[i].getGoodImageSize(scale)
                size = max(tsize, size)

                tsize=psflist[i].getGoodImageSize(scale)
                psf_size = max(tsize, psf_size)

            dims=numpy.array([size,size])
            psf_dims=numpy.array([psf_size,psf_size])


        if self.nrows is not None:
            dims = [self.nrows, self.ncols]

        if self.psf_nrows is not None:
            psf_dims = [self.psf_nrows, self.psf_ncols]

        logger.debug("    image dims: %s" % dims)
        logger.debug("    psf dims:   %s" % psf_dims)

        return dims, psf_dims

    def _get_convolved_object_info(self, objconf):
        """
        get convolved objects and psfs and wcs
        """
        psflist=[]
        cobjlist=[]
        wcslist=[]

        nepoch = objconf.get('nepoch',1)

        object, meta = self._get_object()

        # offset of object, not the epochs
        obj_shift = self._get_obj_shift()
        if obj_shift is not None:
            object = object.shift(
                dx=obj_shift['col_offset'],
                dy=obj_shift['row_offset'],
            )

        if 'randomize_morphology' in objconf:
            raise NotImplementedError("make work with new system "
                                      "of same sized stamps")
            # for multi-epoch sims
            flux = object.getFlux()
            try:
                r50  = object.getHalfLightRadius()
            except:
                r50  = object.calculateHLR()

        for epoch in xrange(nepoch):

            if 'randomize_morphology' in objconf:
                r_flux, r_r50 = self._randomize_morphology(flux, r50)
                object, meta = self._get_object(flux=r_flux, r50=r_r50)

            # this can be random, so only should be called
            # once per object and epoch
            wcs = self._get_galsim_wcs()
            psf, psf_meta = self._get_psf()
            cobj = convolved_object = galsim.Convolve(object, psf)

            wcslist.append(wcs)
            psflist.append(psf)
            cobjlist.append(cobj)

        return cobjlist, psflist, wcslist, meta

    def _do_coadd(self, obslist):
        import psc

        coadd_conf=self['coadd']

        kw={}

        if 'interp' in coadd_conf:
            kw['interp'] = coadd_conf['interp']

        if 'flat_wcs' in coadd_conf:
            kw['flat_wcs'] = coadd_conf['flat_wcs']

        if 'weight_type' in coadd_conf:
            kw['weight_type'] = coadd_conf['weight_type']

        for obs in obslist:
            sigma = numpy.sqrt(1./obs.weight.max())
            obs.noise = self.rng.normal(
                scale=sigma,
                size=obs.image.shape,
            )

        if 'replace_bad_pixels' in coadd_conf:
            # do a fit on the coadd and use it to replace bad pixels
            # in the SE images
            type=coadd_conf['replace_bad_pixels']['type']
            if type == 'me':
                obslist = self._replace_bad_pixels_from_me_fit(obslist)
            elif type=='interp':
                self._replace_bad_pixels_interp(obslist)
            else:
                raise ValueError("unsupported replace_bad_pixels "
                                 "type: '%s'" % type)

        coadder = psc.Coadder(obslist, **kw)

        coadd_obs = coadder.get_coadd()

        coadd_obslist=ngmix.ObsList()
        coadd_obslist.append(coadd_obs)
        coadd_obslist.update_meta_data(obslist.meta)

        if False:
            self._show_coadd(obslist)
            self._show_coadd(coadd_obslist)
        return coadd_obslist

    def _replace_bad_pixels_interp(self, obslist):
        """
        do see bias when noisy
        1 maybe because we need to interpolate the noise also
        2 maybe because the weights were zero and this is messing
          up the noise image in metacal?  No because in coaddsim
          I return a constant weight images

        First trying 1) but also resetting bmask and weight map, so not
        fully controlled.  If it works we can dissect
        """

        coadd_conf=self['coadd']

        assert coadd_conf['use_nsim_noise_image'],"forcing use noise image"
        iconf = coadd_conf['replace_bad_pixels']['interp']
        assert iconf['type']=="cubic","only cubic interpolation for now"

        for obs in obslist:

            im=obs.image
            weight=obs.weight

            if not obs.has_bmask():
                obs.bmask = numpy.zeros(im.shape, dtype='i4')

            bmask = obs.bmask
            noise = obs.noise

            imravel = im.ravel()
            noise_ravel = noise.ravel()
            bmravel = bmask.ravel()
            wtravel = weight.ravel()

            wbad,=numpy.where( (bmravel != 0) | (wtravel == 0.0) )

            if wbad.size > 0:
                #print("        interpolating %d/%d masked or zero weight "
                #      "pixels" % (wbad.size,im.size))

                yy, xx = numpy.mgrid[0:im.shape[0], 0:im.shape[1]]

                x = xx.ravel()
                y = yy.ravel()

                yx = numpy.zeros( (x.size, 2) )
                yx[:,0] = y
                yx[:,1] = x

                wgood, = numpy.where( (bmravel==0) & (wtravel != 0.0) )

                im_interp = self._do_interp(yx, im, wgood, wbad)
                noise_interp = self._do_interp(yx, noise, wgood, wbad)

                obs.image = im_interp
                obs.noise = noise_interp

                # for now set bmask to zero in case downstream is avoiding it
                bmask[:,:]=0

                # maybe want to interpolate weight map too in real data
                obs.weight[:,:] = obs.weight.max()

    def _do_interp(self, yx, im, wgood, wbad):
        import scipy.interpolate

        im_ravel = im.ravel()

        ii = scipy.interpolate.CloughTocher2DInterpolator(
            yx[wgood,:],
            im_ravel[wgood],
            fill_value=0.0,
        )

        im_interp = im.copy()
        im_interp_ravel = im_interp.ravel()

        vals = ii(yx[wbad,:])
        im_interp_ravel[wbad] = vals

        return im_interp


    def _replace_bad_pixels_from_me_fit(self, obslist):
        """
        do a full multi-epoch fit use the model
        to fill in bad pixels in the obslist
        """
        print("replacing bad pixels")
        from ngmix.gexceptions import BootPSFFailure, BootGalFailure
        repconf = self['coadd']['replace_bad_pixels']


        boot = ngmix.bootstrap.Bootstrapper(obslist)
        mconf=repconf['max_pars']
        ppars=repconf['psf_pars']

        scale = obslist[0].jacobian.get_scale()
        psf_Tguess=4.0 * scale**2

        g_prior = ngmix.priors.GPriorBA(0.2, rng=self.rng)
        T_prior= ngmix.priors.FlatPrior(-10.0, 1.e6, rng=self.rng)
        flux_prior= ngmix.priors.FlatPrior(-1.0e+04, 1.0e+09, rng=self.rng)
        cen_prior=ngmix.priors.CenPrior(0.0, 0.0, scale, scale, rng=self.rng)
        prior = ngmix.joint_prior.PriorSimpleSep(
            cen_prior,
            g_prior,
            T_prior,
            flux_prior,
        )

        try:
            boot.fit_psfs(
                ppars['model'],
                psf_Tguess,
                ntry=mconf['ntry'],
                fit_pars=mconf['pars']['lm_pars'],
            )
        except BootPSFFailure:
            raise TryAgainError("failed to fit psf for pixel replace")

        try:

            boot.fit_max(
                repconf['model'],
                mconf['pars'],
                prior=prior,
                ntry=mconf['ntry'],
            )

        except BootGalFailure:
            raise TryAgainError("failed to fit galaxy")

        boot.replace_masked_pixels()

        # get first band
        new_obslist = boot.mb_obs_list[0]
        return new_obslist

    def _add_defects(self, obslist):
        for defect in self['defects']:
            dtype=defect['type']
            if dtype == 'bad_pixel':
                self._add_bad_pixels(obslist, defect)
            elif dtype == 'bad_column':
                self._add_bad_columns(obslist, defect)
            elif dtype=='example_bmasks':
                self._add_example_bmasks(obslist, defect)
            else:
                raise ValueError("bad defect type: '%s'" % defect['type'])

    def _get_example_bmask(self):
        """
        draw randomly from the list of example bmasks

        randomly flip across rows using fliplr and cols using flipud
        """
        nmasks=len(self.bmask_list)
        i = self.rng.randint(0, nmasks)

        bmask = self.bmask_list[i]

        r = self.rng.uniform()
        if r > 0.5:
            bmask = numpy.fliplr(bmask)

        r = self.rng.uniform()
        if r > 0.5:
            bmask = numpy.flipud(bmask)

        return bmask.copy()

    def _add_example_bmasks(self, obslist, defect):
        """
        add single bad pixels with a given rate
        """
        #print("adding example bmasks")
        for obs in obslist:
            if not obs.has_bmask():
                obs.bmask = numpy.zeros( obs.image.shape, dtype='i4' )

            bmask = self._get_example_bmask()

            wbad = numpy.where(bmask != 0)
            if wbad[0].size > 0:
                obs.bmask = bmask
                obs.weight[wbad] = 0.0
                obs.image[wbad] = 1.e9

    def _add_bad_pixels(self, obslist, defect):
        """
        add single bad pixels with a given rate
        """
        print("adding single bad pixels")
        for obs in obslist:
            if not obs.has_bmask():
                obs.bmask = numpy.zeros( obs.image.shape, dtype='i4' )

            bmravel = obs.bmask.ravel()
            imravel = obs.image.ravel()
            wtravel = obs.weight.ravel()

            if defect['rate']=='all':
                num = defect['nper']
            else:
                r = self.rng.uniform()
                if r < defect['rate']:
                    num=1
                else:
                    num=0

            for i in xrange(num):
                # pick a random pixel

                i = self.rng.randint(0, imravel.size)
                #print("    ",i)

                wtravel[i] = 0.0
                imravel[i] = 1.e9
                bmravel[i] = BAD_PIXEL

    def _add_bad_columns(self, obslist, defect):
        """
        add single bad pixels with a given rate
        """
        print("adding bad columns")
        for obs in obslist:
            if not obs.has_bmask():
                obs.bmask = numpy.zeros( obs.image.shape, dtype='i4' )

            bmask  = obs.bmask
            image  = obs.image
            weight = obs.weight

            r = self.rng.uniform()
            if r < defect['rate']:
                # pick a random column

                i = self.rng.randint(0, image.shape[1])
                print("    col:",i)

                weight[:,i] = 0.0
                image[:,i]  = 1.e9
                bmask[:,i]  = BAD_COLUMN

    def _do_coadd_test_straight(self, obslist, type):
        iilist = [] 


        for i,tobs in enumerate(obslist):
            if type=='psf':
                obs = tobs.psf
            else:
                obs = tobs

            jac = obs.jacobian # copy
            wcs = jac.get_galsim_wcs()

            if type=='noise':
                im = obs.noise
            else:
                im = obs.image

            if i==0:
                coadd_im = im.copy()
            else:
                coadd_im += im

        coadd_im /= len(obslist)

        cen = (numpy.array(coadd_im.shape)-1.0)/2.0
        jac.set_cen(row=cen[0], col=cen[1])

        return coadd_im, jac


    def _do_coadd_test_type(self, obslist, type):
        if obslist[0].meta['offset_pixels'] is None:
            return self._do_coadd_test_straight(obslist, type)

        iilist = [] 

        for tobs in obslist:
            if type=='psf':
                #doffset = tobs.meta['psf_offset_pixels']
                doffset = tobs.psf.meta['offset_pixels']
                obs = tobs.psf
                #print("psf offset:",doffset)
            else:
                doffset = tobs.meta['offset_pixels']
                obs = tobs
                #print("obj offset:",doffset)

            if doffset is not None:
                offset=(doffset['col_offset'],doffset['row_offset'])
            else:
                offset=(0.0, 0.0)

            jac = obs.jacobian # copy
            wcs = jac.get_galsim_wcs()

            if type=='noise':
                im = obs.noise
            elif type=='psf':
                im = obs.image
                im = im/im.sum()
            else:
                im = obs.image

            interp=self['coadd']['interp']
            ii = galsim.InterpolatedImage(
                galsim.Image(im, wcs=wcs),
                offset=offset,
                x_interpolant=interp,
            )

            iilist.append(ii)

        # use size and wcs from last one
        ny,nx = im.shape
        coadd_gsim = galsim.Image(nx, ny, wcs=wcs)

        coadd_ii = galsim.Sum(iilist)
        coadd_ii.drawImage(image=coadd_gsim, method='no_pixel')

        coadd_im = coadd_gsim.array
        coadd_im /= len(obslist)

        # assume all same scale
        cen = (numpy.array(obslist[0].image.shape)-1.0)/2.0
        jac.set_cen(row=cen[0], col=cen[1])

        return coadd_im, jac

    def _do_coadd_test(self, obslist):

        coadd_im, jac = self._do_coadd_test_type(obslist, 'image')
        n_im, _ = self._do_coadd_test_type(obslist, 'noise')
        psf_im, psf_jac = self._do_coadd_test_type(obslist, 'psf')

        var = n_im.var()
        weight = n_im*0 + 1.0/var

        psf_obs = ngmix.Observation(
            psf_im,
            weight=obslist[0].psf.weight,
            jacobian=psf_jac,
        )
        coadd_obs = ngmix.Observation(
            coadd_im,
            weight=weight,
            jacobian=jac,
            psf=psf_obs,
        )
        coadd_obs.noise = n_im

        return coadd_obs


    def _show_coadd(self, obslist):
        import biggles
        import images


        for i,obs in enumerate(obslist):
            if len(obslist) > 1:
                title='image %d' % i
            else:
                title='coadd'

            tab = biggles.Table(2, 1)
            ppsf = images.multiview(obs.psf.image, title='psf',width=1000,height=1000,show=False)
            pim = images.multiview(obs.image, title=title,width=1000,height=1000,show=False)
            tab[0,0] = ppsf
            tab[1,0] = pim
            tab.show()

        if raw_input('hit a key (q to quit): ')=='q':
            stop

    def _randomize_morphology(self, flux, r50):

        rng=self.rng

        rspec = self['object']['randomize_morphology']

        if rspec['type']=='uniform':
            low,high=rspec['range']
            r_flux = flux *(1.0 + rng.uniform(low=low, high=high))
            r_r50  = r50  *(1.0 + rng.uniform(low=low, high=high))
        else:
            raise ValueError("randomize_morphology must be of "
                             "type uniform, got '%s'" % rspec['type'])

        return r_flux, r_r50

    def _get_obs(self, psf, cobj, wcs, dims, psf_dims):

        offset_pixels = self._get_epoch_offset()

        noise_obj=self._get_noise()

        if self['psf']['shift_psf']:
            psf_offset_pixels = offset_pixels
        else:
            psf_offset_pixels = None

        #print("offset:",offset_pixels)
        #print("offset psf:",psf_offset_pixels)
        psf_im, psf_jacob = self._get_psf_image(
            psf,
            wcs,
            psf_dims,
            noise_obj,
            psf_offset_pixels,
        )

        obj_im, obj_im_orig, obj_jacob =  self._get_object_image(
            cobj,
            wcs,
            dims,
            noise_obj,
            offset_pixels,
        )

        ivar = 1.0/noise_obj.sigma**2
        psf_weight = numpy.zeros( psf_im.shape ) + ivar
        obj_weight = numpy.zeros( obj_im.shape ) + ivar

        psf_obs = ngmix.observation.Observation(
            psf_im,
            weight=psf_weight,
            jacobian=psf_jacob,
        )

        obs = ngmix.observation.Observation(
            obj_im,
            weight=obj_weight,
            jacobian=obj_jacob,
            psf=psf_obs,
        )
        obs.image_orig = obj_im_orig
        obs.update_meta_data({'offset_pixels':offset_pixels})
        obs.psf.update_meta_data({'offset_pixels':psf_offset_pixels})

        return obs

    def _make_biased_psf(self,psf,psf_dims):
        """
        produce a biased psf

        currently support dilation
        """
        pb = self['psf_bias']
        if 'dilate' in pb:
            dilate = pb['dilate']
            #print("dilating psf by:",dilate)
            psf = psf.dilate(dilate)
            psf_dims = ( psf_dims * dilate ).astype('i4')
        else:
            raise RuntimeError("expected psf bias parameters")

        return psf, psf_dims

    def _get_psf_image(self, psf, wcs, psf_dims, noise_obj, offset):
        """
        """

        if 'psf_bias' in self:
            psf, psf_dims = self._make_biased_psf(psf,psf_dims)

        gsimage = self._make_gsimage(
            psf,
            wcs,
            nrows=psf_dims[0],
            ncols=psf_dims[1],
            offset=offset,
        )

        gsimage.addNoiseSNR(
            noise_obj,
            self['psf']['s2n'],
            preserve_flux=True,
        )


        image = gsimage.array

        # no offset, jacobian is straightforward
        dims = numpy.array(gsimage.array.shape)
        row, col = (numpy.array(dims)-1.0)/2.0
        if offset is not None:
            row += offset['row_offset']
            col += offset['col_offset']

        jacob = self._get_jacobian(wcs, row, col)

        return image, jacob


    def _get_object_image(self, convolved_object, wcs, dims, noise_obj, offset):
        """
        convolve the 
        """

        gsimage = self._make_gsimage(
            convolved_object,
            wcs,
            nrows=dims[0],
            ncols=dims[1],
            offset=offset,
        )


        # find centroid and get the jacobian
        # for centroid finding, use original image before
        # noise

        image_orig = gsimage.array.copy()

        if self['use_canonical_center']:
            tdims = numpy.array(image_orig.shape)
            row, col = (tdims-1)/2.0

            if offset is not None:
                row += offset['row_offset']
                col += offset['col_offset']
            logger.debug("using canonical center %s" % (row,col))
        else:
            row, col = find_centroid(image_orig, self.rng, offset=offset)

        jacob = self._get_jacobian(wcs, row, col)

        # add noise
        gsimage.addNoise(noise_obj)

        image = gsimage.array

        return image, image_orig, jacob

    def _get_jacobian(self, wcs, row, col):
        """
        use converter from galsim wcs to jacobian
        """
        jacob = ngmix.Jacobian(
            wcs=wcs,
            row=row,
            col=col,
        )

        return jacob

    def _make_gsimage(self,
                      gs_obj,
                      wcs,
                      nrows=None,
                      ncols=None,
                      offset=None):
        """
        if nrows,ncols None, dims are chosen by galsim
        """

        gs_offset=offset
        if offset is not None:
            gs_offset = (offset['col_offset'], offset['row_offset'])

        return gs_obj.drawImage(
            nx=ncols,
            ny=nrows,
            wcs=wcs,
            dtype=numpy.float64,
            offset=gs_offset,
        )


    def _get_galsim_wcs(self):
        """
        set basic wcs info

        The actual ngmix jacobian will be created from this later
        """
        from math import cos, sin

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
            elif 'rotate' in ws:
                if ws['rotate']['type'] == 'uniform':
                    jmatrix = numpy.array(
                        [[dudx,dudy],
                         [dvdx,dvdy]]
                    )
                    angle = self.rng.uniform(
                        low=0.0,
                        high=numpy.pi*2,
                    )
                    angle=1.0*numpy.pi/1.0
                    rotmatrix = numpy.array(
                        [[ cos(angle), -sin(angle) ],
                         [ sin(angle),  cos(angle) ]]
                    )

                    rotated_jmatrix = numpy.dot(rotmatrix,jmatrix)
                    dudx = rotated_jmatrix[0,0]
                    dudy = rotated_jmatrix[0,1]
                    dvdx = rotated_jmatrix[1,0]
                    dvdy = rotated_jmatrix[1,1]

                    """
                    import images
                    print("original")
                    images.imprint(jmatrix)
                    print("rotated",angle*180.0/numpy.pi)
                    images.imprint(rotated_jmatrix)
                    print("-"*70)
                    stop
                    """
        else:
            dudx=1.0
            dudy=0.0
            dvdx=0.0
            dvdy=1.0

        return galsim.JacobianWCS(
            dudx,
            dudy,
            dvdx,
            dvdy,
        )

    def _get_object(self, **kw):
        """
        get the galsim representation of the object

        parameters
        ----------
        flux: float, optional
            Force the given flux
        r50: float, optional
            Force the given r50

        returns
        -------
        (gsobj, meta)
           The galsim objects and metadata in a dictionary
        """

        obj, meta  = self.object_maker(**kw)

        if self.shear_pdf is not None:
            shear, shindex = self.shear_pdf.get_shear()
            obj = obj.shear(g1=shear.g1, g2=shear.g2)

            meta['shear'] = (shear.g1, shear.g2)
            meta['shear_index'] = shindex

        return obj, meta

    def _get_psf(self):
        """
        get the galsim representation of the psf
        """
        return self.psf_maker()

    def _set_pdfs(self):

        self._set_noise()
        self._set_sizes()
        self._set_obj_shift()
        self._set_epoch_offset()

    def _set_sizes(self):
        if 'stamp_size' in self['object']:
            self.nrows,self.ncols=self['object']['stamp_size']
        else:
            self.nrows,self.ncols=None,None

        if 'stamp_size' in self['psf']:
            self.psf_nrows,self.psf_ncols=self['psf']['stamp_size']
        else:
            self.psf_nrows,self.psf_ncols=None,None

    def _get_noise(self):
        """
        choose one of the noises at random
        """
        return random.choice(self.noise_objects)

    def _set_noise(self):
        if 'noise' not in self:
            raise ValueError("set noise in the ['images'] "
                             "section of the config")

        if 's2n' not in self['psf']:
            raise ValueError("set psf s2n in the ['images']['psf'] "
                             "section of the config")

        if isinstance(self['noise'], dict):
            assert self['noise']['type']=='uniform'
            nrand = self['noise']['num']
            sigma_min,sigma_max = self['noise']['sigma_range']

            sigmas=self.rng.uniform(
                low=sigma_min,
                high=sigma_max,
                size=nrand,
            )


        else:
            sigmas = [self['noise']]

        self.noise_objects = [self._make_one_noise(s) for s in sigmas]

    def _make_one_noise(self, sigma):
        return galsim.GaussianNoise(
            self.galsim_rng,
            sigma,
        )


    def _get_obj_shift(self):
        """
        shift of object on the sky
        """
        return self._get_offset_generic(self.obj_shift_pdf)

    def _get_epoch_offset(self):
        """
        offset of an epoch, not object
        """
        return self._get_offset_generic(self.epoch_offset_pdf)

    def _get_offset_generic(self, pdf):
        """
        the shift of the object, not the epoch shift
        """
        if pdf is not None:
            if hasattr(pdf,'sample2d'):
                coff1,coff2 = pdf.sample2d()
            else:
                coff1 = pdf.sample()
                coff2 = pdf.sample()

            offset=(coff1,coff2)
            offset = {
                'row_offset':coff1,
                'col_offset':coff2,
            }
        else:
            offset=None

        return offset


    def _set_obj_shift(self):
        cr=self['object'].get('cen_shift',None)
        self.obj_shift_pdf = self._get_shift_pdf(cr)

    def _set_epoch_offset(self):
        cr=self['offset']
        self.epoch_offset_pdf = self._get_shift_pdf(cr)

    def _get_shift_pdf(self, cr):
        if cr is None:
            pdf=None
        else:
            type=cr.get('type','uniform')
            if type=='uniform':
                pdf=ngmix.priors.FlatPrior(
                    -cr['radius'], cr['radius'],
                    rng=self.rng,
                )
            elif type=='disk':
                pdf=ngmix.priors.ZDisk2D(
                    cr['radius'],
                    rng=self.rng,
                )

            elif type=='annulus':
                pdf=ngmix.priors.ZAnnulus(
                    cr['rmin'],
                    cr['rmax'],
                    rng=self.rng,
                )


            else:
                raise ValueError("cen shift type should be 'uniform'")

        return pdf

    def _load_bmasks(self):
        """
        load masks from a multi-extension fits file

        if add_rotated is set, a 90 degree rotated version
        is added to cancel symmetries in the mask, such as
        bad columns
        """

        mask_file=None
        if 'defects' in self:
            for defect in self['defects']:
                if defect['type']=='example_bmasks':
                    mask_file=os.path.expandvars(defect['file'])
                    add_rotated=defect.get('add_rotated',False)
                    break

        if mask_file is None:
            print("no bmasks to load")
            return

        print("Loading masks from:",mask_file)

        bmask_list=[]
        with fitsio.FITS(mask_file) as fits:

            for hdu in fits:
                if hdu.get_extname()=='catalog':
                    continue

                mask = hdu.read()
                if add_rotated:
                    #print("adding rotated")
                    rm = numpy.rot90(mask)
                    mask = mask + rm

                bmask_list.append( mask )

        print("    loaded %d masks" % len(bmask_list))
        self.bmask_list=bmask_list


class NbrObservationMaker(ObservationMaker):
    """

    In this version we generate different objects totally separately and then
    just add the images. I think this might not be right because then the whole
    thing might not be sheared in a self-consistent way

    Get an obs with neighbors

    returns
    -------
    ngmix.Observation:

    notes
    ------
        - the noise gets increased since the images are added
        - you should fix the stamp size so the images can be added
        - the cen_shift does not apply to the first object drawn
          only subsequent ones

        - this requires the additional config parameter
            object:
                nnbr:   2 # number of neighbors
                nepoch: 1 # only one epoch for now
        - number of epochs must be 1 for now
    """
    def __call__(self):


        objconf=self['object']
        assert objconf['nepoch']==1

        nobject = 1 + objconf['nnbr']
        for i in xrange(nobject):
            logger.debug("    object %d of %d" % (i+1,nobject))
            if i==0:
                self._save_cen_pdf()

            obslist = super(NbrObservationMaker,self).__call__()
            obs=obslist[0]

            if i==0:
                self._restore_cen_pdf()

                new_obslist = deepcopy(obslist)


                image  = obs.image
                weight = obs.weight
                var    = 1.0/weight
            else:
                image += obs.image
                var   += 1.0/obs.weight

        weight = 1.0/var
        obs = new_obslist[0]

        obs.image_orig = obs.image
        obs.weight_orig = obs.weight

        obs.image = image.copy()
        obs.weight = weight.copy()

        if False:
            import images
            #images.view(obs.image, width=800,height=800)
            images.multiview(obs.image, width=800,height=800)
            if 'q'==raw_input('hit a key: '):
                stop
            # just to keep things going
            #raise TryAgainError("for testing")

        return new_obslist

    def _restore_cen_pdf(self):
        self.obj_shift_pdf = self.obj_shift_pdf_saved

    def _save_cen_pdf(self):
        self.obj_shift_pdf_saved = self.obj_shift_pdf
        self.obj_shift_pdf = None


class NbrObservationMakerMulti(ObservationMaker):
    """
    multi means an observation is returned with center
    on each of the objects
    """
    def __init__(self, *args, **kw):
        super(NbrObservationMakerMulti,self).__init__(*args, **kw)

        self._set_nbr_sky_shift()
        self._set_nbr_size_dilation()

    def __call__(self):
        
        obslist=super(NbrObservationMakerMulti,self).__call__()
        
        obs = obslist[0]
        jac = obs.jacobian
        meta=obslist.meta

        shiftlist=meta['shiftlist']
        nobj = len(shiftlist)
        logger.debug("nobj: %d" % nobj)

        allobs = [] 

        for i in xrange(nobj):


            if i==0:
                allobs.append( obslist )
            else:
                # this works because we demand trivial wcs
                dskyrow, dskycol = shiftlist[i]

                tjac = jac.copy()
                scale = jac.get_scale()
                row,col = jac.get_cen()
                drow,dcol = dskyrow/scale, dskycol/scale

                row0 = row + drow
                col0 = col + dcol
                logger.debug("new cen:" % (row0,col0))
                tjac.set_cen(row=row0, col=col0)

                tobs = ngmix.Observation(
                    obs.image.copy(),
                    obs.weight.copy(),
                    jacobian=tjac,
                    meta=obs.meta,
                    psf=obs.psf,
                )

                tlist = ngmix.ObsList()
                tlist.append(tobs)
                tlist.update_meta_data(obslist.meta)

                allobs.append(tlist)


        if False:
            import images
            #images.view(obs.image, width=800,height=800)
            obs = allobs[0][0]
            images.multiview(obs.image, width=800,height=800)
            if 'q'==raw_input('hit a key: '):
                stop
            # just to keep things going
            #raise TryAgainError("for testing")
        return allobs



    def _get_object(self, **kw):
        """
        get the galsim representation of the object

        parameters
        ----------
        flux: float, optional
            Force the given flux
        r50: float, optional
            Force the given r50

        returns
        -------
        (gsobj, meta)
           The galsim objects and metadata in a dictionary
        """

        nnbr = self['nnbr']

        objlist=[]
        shiftlist = []
        parlist=[]

        nobj = 1 + nnbr
        for i in xrange(nobj):
            tobj, tmeta  = self.object_maker(**kw)

            shift=None
            if i > 0:

                # always dilate first
                dilation = self._get_nbr_dilation()
                if dilation is not None:
                    tobj = tobj.dilate(dilation)
                    logger.debug("dilation: %s flux: %s" % (dilation,tobj.getFlux()))

                if self.nbr_sky_shift_pdf is not None:
                    shift = self._get_nbr_sky_shift()
                    tobj = tobj.shift(dx=shift[1], dy=shift[0])
            else:
                meta=tmeta

            objlist.append(tobj)
            shiftlist.append(shift)

            pars={}
            pars.update(tmeta)
            parlist.append(pars)

        obj = galsim.Add(objlist)

        meta['shiftlist'] = shiftlist
        meta['parlist'] = parlist

        if self.shear_pdf is not None:
            shear, shindex = self.shear_pdf.get_shear()
            obj = obj.shear(g1=shear.g1, g2=shear.g2)

            meta['shear'] = (shear.g1, shear.g2)
            meta['shear_index'] = shindex

        return obj, meta

    def _get_nbr_sky_shift(self):
        sky_shift_pdf = self.nbr_sky_shift_pdf
        if hasattr(sky_shift_pdf,'sample2d'):
            coff1,coff2 = sky_shift_pdf.sample2d()
        else:
            coff1 = sky_shift_pdf.sample()
            coff2 = sky_shift_pdf.sample()

        offset=(coff1,coff2)

        logger.debug("sky shift: %g,%g" % offset)

        return offset

    def _get_galsim_wcs(self):
        """
        set basic wcs info

        The actual ngmix jacobian will be created from this later
        """
        if 'wcs' in self:
            raise ValueError("only trivial wcs for nbrs")

        dudx=1.0
        dudy=0.0
        dvdx=0.0
        dvdy=1.0

        return galsim.JacobianWCS(
            dudx,
            dudy,
            dvdx,
            dvdy,
        )

    def _get_nbr_dilation(self):
        return self._dilation

    def _set_nbr_size_dilation(self):
        spec = self['nbr_size_ratio']
        if spec is not None:
            self._dilation = spec
        else:
            raise  ValueError("currently only support "
                              "constant for size ratio")

    def _set_nbr_sky_shift(self):
        cr=self.get('nbr_sky_shift',None)

        if cr is None:
            self.nbr_sky_shift_pdf=None
        else:
            type=cr.get('type','uniform')
            if type=='uniform':
                self.nbr_sky_shift_pdf=ngmix.priors.FlatPrior(
                    -cr['radius'], cr['radius'],
                    rng=self.rng,
                )
            elif type=='disk':
                self.nbr_sky_shift_pdf=ngmix.priors.ZDisk2D(
                    cr['radius'],
                    rng=self.rng,
                )

            elif type=='annulus':
                self.nbr_sky_shift_pdf=ngmix.priors.ZAnnulus(
                    cr['rmin'],
                    cr['rmax'],
                    rng=self.rng,
                )


            else:
                raise ValueError("cen shift type should be 'uniform'")




def find_centroid(image, rng, offset=None, maxiter=200, ntry=4):
    """
    use AM to fit a single gaussian
    """


    dims = numpy.array(image.shape)
    row0, col0 = (dims-1)/2.0

    if offset is not None:
        row0 += offset['row_offset']
        col0 += offset['col_offset']


    if False:
        import images
        images.multiview(image, width=800,height=800)
        if 'q'==raw_input('hit a key: '):
            stop



    j=ngmix.UnitJacobian(row=row0, col=col0)

    obs = ngmix.Observation(image, jacobian=j)

    guess_T = 4.0
    for i in xrange(ntry):
        fitter = ngmix.admom.run_admom(obs, guess_T, rng=rng)

        res=fitter.get_result()
        if res['flags']==0:
            break

    if res['flags'] != 0:
        if False:
            import images
            images.view(image, width=800,height=800)
            if 'q'==raw_input('hit a key: '):
                stop

        raise TryAgainError("could not fit 1 gauss to get centroid")

    pars=res['pars']
    row=pars[0]
    col=pars[1]

    row = row0 + row
    col = col0 + col

    return row,col

def get_expected_s2n(obslist):
    """
    maximal s2n
    """

    sum = 0.0

    for obs in obslist:
        sum += (obs.image_orig**2 *obs.weight).sum()

    s2n = numpy.sqrt(sum)
    return s2n


