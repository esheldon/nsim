from __future__ import print_function
try:
    xrange
except:
    xrange=range

import numpy
import galsim
import ngmix
from . import pdfs

from .util import TryAgainError

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
        nepoch: 10

        stamp_size: [48,48] # optional
        cen_shift:
            type: "uniform"
            radius: 0.5 # pixels
    """
    def __init__(self, config, psf_maker, object_maker, rng):
        self.update(config)

        self.psf_maker=psf_maker
        self.object_maker=object_maker
        self.rng=rng

        self._set_pdfs()

    def __call__(self):

        objconf=self['object']

        object, meta = self._get_object()


        if 'randomize_morphology' in objconf:
            flux = object.getFlux()
            try:
                r50  = object.getHalfLightRadius()
            except:
                r50  = object.calculateHLR()

        obslist = ngmix.observation.ObsList()
        for epoch in xrange(objconf['nepoch']):

            psf, psf_meta = self._get_psf()

            if 'randomize_morphology' in objconf and epoch > 0:
                # random ellipticity but same flux and size
                r_flux, r_r50 = self._randomize_morphology(flux, r50)
                object, meta = self._get_object(flux=r_flux, r50=r_r50)

            obs = self._get_obs(psf, object)

            obslist.append( obs )

        meta['s2n']  = get_expected_s2n(obslist)

        obslist.update_meta_data(meta)
        return obslist

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

    def _get_obs(self, psf, object):

        wcs=self._get_galsim_wcs()

        psf_im, psf_jacob = self._get_psf_image(
            psf,
            wcs,
        )

        convolved_object = galsim.Convolve(object, psf)
        obj_im, obj_im_orig, obj_jacob = self._get_object_image(
            convolved_object,
            wcs,
        )

        ivar = 1.0/self['noise']**2
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

        return obs


    def _get_psf_image(self, psf, wcs):
        """
        """
        gsimage = self._make_gsimage(
            psf,
            wcs,
            nrows=self.psf_nrows,
            ncols=self.psf_ncols,
        )

        self._add_psf_noise(gsimage)
        image = gsimage.array

        # no offset, jacobian is straightforward
        dims = numpy.array(gsimage.array.shape)
        row, col = (numpy.array(dims)-1.0)/2.0

        jacob = self._get_jacobian(wcs, row, col)

        return image, jacob


    def _get_object_image(self, convolved_object, wcs):
        """
        convolve the 
        """

        offset=self._get_offset()
        gsimage = self._make_gsimage(
            convolved_object,
            wcs,
            nrows=self.nrows,
            ncols=self.ncols,
            offset=offset,
        )

        image_orig = gsimage.array.copy()

        self._add_object_noise(gsimage)
        image=gsimage.array

        row, col = find_centroid(image, self.rng)
        jacob = self._get_jacobian(wcs, row, col)

        return image, image_orig, jacob


    def _add_psf_noise(self,image):
        image.addNoiseSNR(
            self.gaussian_image_noise,
            self['psf']['s2n'],
        )

    def _add_object_noise(self,image):
        image.addNoise(
            self.gaussian_image_noise,
        )


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
        return gs_obj.drawImage(
            nx=ncols,
            ny=nrows,
            wcs=wcs,
            dtype=numpy.float64,
            offset=offset,
        )


    def _get_galsim_wcs(self):
        """
        set basic wcs info

        The actual ngmix jacobian will be created from this later
        """
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
        """
        return self.object_maker(**kw)

    def _get_psf(self):
        """
        get the galsim representation of the psf
        """
        return self.psf_maker()

    def _set_pdfs(self):

        self._set_noise()
        self._set_sizes()
        self._set_offsets()

    def _set_sizes(self):
        if 'stamp_size' in self['object']:
            self.nrows,self.ncols=self['object']['stamp_size']
        else:
            self.nrows,self.ncols=None,None

        if 'stamp_size' in self['psf']:
            self.psf_nrows,self.psf_ncols=self['psf']['stamp_size']
        else:
            self.psf_nrows,self.psf_ncols=None,None


    def _set_noise(self):
        if 'noise' not in self:
            raise ValueError("set noise in the ['images'] "
                             "section of the config")

        if 's2n' not in self['psf']:
            raise ValueError("set psf s2n in the ['images']['psf'] "
                             "section of the config")

        self.galsim_rng = galsim.BaseDeviate(self.rng.randint(0,2**30))
        self.gaussian_image_noise=galsim.GaussianNoise(
            self.galsim_rng,
            self['noise'],
        )

    def _get_offset(self):
        if self.cen_pdf is not None:
            coff1 = self.cen_pdf.sample()
            coff2 = self.cen_pdf.sample()

            offset=(coff1,coff2)
        else:
            offset=None

        if offset is not None:
            print("    offset: %g,%g" % offset)

        return offset


    def _set_offsets(self):
        cr=self['object'].get('cen_shift',None)

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


def find_centroid(image, rng, maxiter=200, ntry=4):
    """
    use AM to fit a single gaussian
    """

    dims = numpy.array(image.shape)
    row0, col0 = (dims-1)/2.0

    j=ngmix.UnitJacobian(row=row0, col=col0)

    obs = ngmix.Observation(image, jacobian=j)

    guess_T = 4.0
    for i in xrange(ntry):
        fitter = ngmix.admom.run_admom(obs, guess_T, rng=rng)

        res=fitter.get_result()
        if res['flags']==0:
            break

    if res['flags'] != 0:
        raise TryAgainError("could not fit 1 gauss")

    pars=res['pars']
    row=pars[0]
    col=pars[1]

    row = row0 + row
    col = col0 + col
    print("    center:",row,col)

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


