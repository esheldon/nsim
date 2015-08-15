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
        base_rng = galsim.BaseDeviate(seed)

        self.gauss_noise = galsim.GaussianNoise(base_rng, sigma=self['skysig'])

        self.base_rng = base_rng

        pprint.pprint(self)

    def get_image_pair(self):
        """
        Get an image pair, with noise added
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

        gal, psf = self.get_galsim_objects()

        psf_obs = self._make_obs(psf, self['psf']['s2n'])
        gal_obs = self._make_obs(gal, s2n)

        gal_obs.set_psf(psf_obs)

        return {'obs':gal_obs}

    def _make_obs(self, gs_obj, s2n):

        dims = self['stamp_size']
        gsimage = galsim.ImageD(dims[0], dims[1])

        gs_obj.drawImage(gsimage, scale=self['pixel_scale'])

        image_ref = gsimage.array

        jacob = self.get_jacobian(image_ref)

        gsimage.addNoiseSNR(self.gauss_noise, s2n)
        
        image = gsimage.array.copy()

        weight = numpy.zeros( image.shape ) 1.0/self['skysig']**2

        obs = ngmix.Observation(image, weight=weight, jacobian=jacob)
        return obs

    def get_jacobian(self, image):
        """
        find the best center and set the jacobian center there
        """
        fitter = quick_fit_gauss(image)
        row,col = fitter.get_gmix().get_cen()

        scale=self['pixel_scale']
        return ngmix.Jacobian(row, col,
                              scale,
                              0.0,
                              0.0,
                              scale)

    def get_galsim_objects(self):
        """
        get the galaxy and psf galsim objects
        """

        coff1 = self.cen_pdf.sample()
        coff2 = self.cen_pdf.sample()

        cenoff=(coff1,coff2)

        psf  = self.get_psf_obj(cenoff)
        gal0 = self.get_gal_obj(cenoff)

        gal = galsim.Convolve([psf, gal0])

        return gal, psf

    def get_gal_obj(self, cenoff):
        """
        get the galsim object for the galaxy model
        """

        pars=self.get_galaxy_pars()

        r50 = pars['r50']
        bd_ratio = pars['bd_ratio']
        s1,s2=self['shear']

        disk_flux = bd_ratio-1.0
        bulge_flux = bd_ratio

        disk = galsim.Exponential(flux=disk_flux, half_light_radius=r50)
        bulge = galsim.DeVaucouleurs(flux=bulge_flux, half_light_radius=r50)

        disk = disk.shear(g1=s1, g2=s2)
        bulge = bulge.shear(g1=s1, g2=s2)

        disk = disk.shift(dx=cenoff[0], dy=cenoff[1])

        boff = pars['boff']
        bulge = bulge.shift(dx=boff[0], dy=boff[1])

        gal = galsim.Add([disk, bulge])


        return gal

    def get_psf_obj(self, cenoff):
        """
        get the galsim object for the psf
        """

        pspec = self['psf']

        psf = galsim.Moffat(beta=pspec['beta'],
                            half_light_radius=pspec['r50'])
        psf = psf.shear(g1=pspec['shape'][0],
                        g2=pspec['shape'][1])

        psf = psf.shift(dx=cenoff[0], dy=cenoff[1])
        return psf

    def get_galaxy_pars(self):
        """
        Get pair parameters

        if not random, then the mean pars are used, except for cen and g1,g2
        which are zero
        """


        g1,g2 = self.g_pdf.sample2d()
        g1=g1[0]
        g2=g2[0]

        r50 = self.r50_pdf.sample()

        bd_ratio = self.bdratio_pdf.sample()

        # distribution is in units of r50
        boff = r50*self.bulge_offset_pdf.sample2d()

        return {'g':(g1,g2),
                'r50':r50,
                'bd_ratio':bd_ratio,
                'boff':boff}

    def _set_pdfs(self):
        """
        Set all the priors
        """
        import ngmix

        self['pixel_scale'] = self.get('pixel_scale',1.0)

        omodel=self['obj_model']

        s2n_r = omodel['s2n']['range']
        self.s2n_pdf=ngmix.priors.FlatPrior(s2n_r[0], s2n_r[1])

        g_spec=omodel['g']
        self.g_pdf=ngmix.priors.GPriorBA(g_spec['sigma'])

        r50_r = omodel['r50']['range']
        self.r50_pdf=ngmix.priors.FlatPrior(r50_r[0], r50_r[1])

        cr=omodel['cen_shift']['range']
        self.cen_pdf=ngmix.priors.FlatPrior(cr[0], cr[1])

        bs_spec=omodel['bulge_shift']
        self.bulge_offset_pdf = ngmix.priors.ZDisk2D(bs_spec['radius'])

        bdr = omodel['bd_ratio']['range']
        self.bdratio_pdf=ngmix.priors.FlatPrior(bdr[0], bdr[1])


    
def quick_fit_gauss(image, maxiter=4000, tol=1.0e-6, ntry=4):
    """
    use EM to fit a single gaussian
    """

    dims = numpy.array(image.shape)
    cenguess = (dims-1)/2.0

    j=ngmix.UnitJacobian(cenguess[0], cenguess[1])

    obs = ngmix.Observation(image, jacobian=j)

    guess_T = 4.0

    for i in xrange(ntry):
        guess_pars = [0.0 + 0.1*srandu(),
                      0.0 + 0.1*srandu(),
                      0.0 + 0.02*srandu(),
                      0.0 + 0.02*srandu(),
                      guess_T*(1.0 + 0.05*srandu()),
                      1.0 + 0.05*srandu() ]
        fitter=ngmix.em.fit_em(obs, guess, maxiter=maxiter, tol=tol)

        res=fitter.get_result()
        if res['flags']==0:
            break
    
    if res['flags'] != 0:
        raise TryAgainError("could not fit 1 gauss")

    return fitter
