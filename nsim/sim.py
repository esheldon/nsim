"""
Simulate images and fit them.

Currently importing ngmix locally because it is so slow; numba does not yet
cache object files. Will make global import when caching is part of numba

additional dependence on
    emcee for MCMC fitting and
    fitsio if checkpointing

for example sim and run configs, see the config subdirectory

"""
from __future__ import print_function
import os
import time
import pprint

import numpy
from numpy import array, zeros, ones, log, log10, exp, sqrt, diag
from numpy.random import random as randu
from numpy.random import randn

import ngmix
from ngmix.fitting import print_pars
from ngmix.observation import Observation
from ngmix.gexceptions import GMixRangeError 

# image will be 2*5 sigma on a side
NSIGMA_IMAGE=5.0

class NGMixSim(dict):
    def __init__(self, sim_conf, s2n):
        """
        Simulate images

        example
        -------
        s2n=1000.0
        sim=NGMixSim(conf, s2n)

        for i in xrange(1000):
            impair = si.get_image_pair()
            # process images
        """

        # seed a new MT random generator from devrand
        # and return the object

        seed_global_devrand()

        self._set_config(sim_conf)

        # this might be mean for correspond to a flux threshold, 
        # depending on the type of run
        self['s2n']=s2n

        self._set_shear()

        self._create_model_set()

        self._set_pdf()
        self._make_psf()

        self._set_noise()

        pprint.pprint(self)

    def get_image_pair(self):
        """
        Get an image pair, with noise added
        """

        imdict=self.get_image_pair_nonoise()

        obs1_0 = imdict['im1']['obs0']
        obs2_0 = imdict['im2']['obs0']

        im1_0 = obs1_0.image
        im2_0 = obs2_0.image
        im1,nim1=self._add_noise(im1_0)
        im2,nim2=self._add_noise(im2_0)

        wt=numpy.zeros(im1_0.shape) + self['ivar']

        obs1 = Observation(im1,
                           weight=wt,
                           psf=obs1_0.psf,
                           jacobian=obs1_0.jacobian)
        obs2 = Observation(im2,
                           weight=wt,
                           psf=obs2_0.psf,
                           jacobian=obs2_0.jacobian)

        imdict['im1']['obs'] = obs1
        imdict['im2']['obs'] = obs2
        imdict['im1']['nim'] = nim1
        imdict['im2']['nim'] = nim2

        im1_s2n = self._get_model_s2n(im1_0)
        im2_s2n = self._get_model_s2n(im2_0)

        imdict['im1']['s2n'] = im1_s2n
        imdict['im2']['s2n'] = im2_s2n

        imdict['im1']['noise'] = self['skysig']
        imdict['im2']['noise'] = self['skysig']

        return imdict



    def get_image_pair_nonoise(self, random=True):
        """
        get a model image

        If random is True, use draw random values from the priors.
        Otherwise use the mean of the priors

        """
        import ngmix

        while True:
            try:
                pars1, pars2, p1_noshear, p2_noshear = self.get_pair_pars(random=random)

                obj_model = self.model_sampler()
                print("sim model:",obj_model)
                #print_pars(pars1,front="pars1:")
                #print_pars(pars2,front="pars2:")
                #print("psf T:",self.psf_gmix_true.get_T())

                gm1_pre=ngmix.gmix.GMixModel(pars1, obj_model)
                gm2_pre=ngmix.gmix.GMixModel(pars2, obj_model)

                gm1  = gm1_pre.convolve(self.psf_gmix_true)
                gm2  = gm2_pre.convolve(self.psf_gmix_true)

                # in case not doing a ring
                T = max(gm1.get_T(), gm2.get_T())

                dims_pix, cen0_pix = self._get_dims_cen(T + self['psf_T'])

                # conversion between pixels and sky in arcsec
                self.jacobian=ngmix.jacobian.Jacobian(cen0_pix[0],
                                                      cen0_pix[1],
                                                      self['pixel_scale'],
                                                      0.0,
                                                      0.0,
                                                      self['pixel_scale'])

                psf_cen=pars1[0:0+2].copy()
                #psf_cen=[0,0]
                psf_image=self._get_psf_image(dims_pix, psf_cen)

                psf_obs = Observation(psf_image, jacobian=self.jacobian)

                nsub = self['nsub']
                im1=gm1.make_image(dims_pix, nsub=nsub, jacobian=self.jacobian)
                im2=gm2.make_image(dims_pix, nsub=nsub, jacobian=self.jacobian)

                obs1 = Observation(im1, psf=psf_obs, jacobian=self.jacobian)
                obs2 = Observation(im2, psf=psf_obs, jacobian=self.jacobian)
                
                out={'im1':{'pars':pars1,
                            'pars_noshear':p1_noshear,
                            'gm_pre':gm1_pre,
                            'gm':gm1,
                            'obs0':obs1,
                            'model':obj_model},
                     'im2':{'pars':pars2,
                            'pars_noshear':p2_noshear,
                            'gm_pre':gm2_pre,
                            'gm':gm2,
                            'obs0':obs2,
                            'model':obj_model},
                    }

                break
            except GMixRangeError as err:
                print("caught range exception: %s" % str(err))

        return out

    def get_pair_pars(self, random=True):
        """
        Get pair parameters

        if not random, then the mean pars are used, except for cen and g1,g2
        which are zero
        """
        import ngmix
        from numpy import pi

        if random:
            pars1 = self.pdf.sample()

            if self['do_ring']:
                # use same everything but rotated 90 degrees
                pars2=pars1.copy()

                g1=pars1[2]
                g2=pars1[3]
            
                shape2=ngmix.shape.Shape(g1,g2)
                shape2.rotate(pi/2.0)

                pars2[2]=shape2.g1
                pars2[3]=shape2.g2
            else:
                # use different ellipticity
                pars2 = self.pdf.sample()

        else:
            samples=self.pdf.sample(10000)
            pars1 = samples.mean(axis=0)

            pars1[0:0+4] = 0.0
            pars2=pars1.copy()

        pars1_noshear=pars1.copy()
        pars2_noshear=pars2.copy()

        shape1=ngmix.shape.Shape(pars1[2],pars1[3])
        shape2=ngmix.shape.Shape(pars2[2],pars2[3])

        Tround1 = ngmix.moments.get_Tround(pars1[4], pars1[2], pars1[3])
        Tround2 = ngmix.moments.get_Tround(pars2[4], pars2[2], pars2[3])

        shear=self._get_shear()
        shape1.shear(shear[0], shear[1])
        shape2.shear(shear[0], shear[1])

        Tsheared1 = ngmix.moments.get_T(Tround1, shape1.g1, shape1.g2)
        Tsheared2 = ngmix.moments.get_T(Tround2, shape2.g1, shape2.g2)

        pars1[2]=shape1.g1
        pars1[3]=shape1.g2
        pars2[2]=shape2.g1
        pars2[3]=shape2.g2

        pars1[4] = Tsheared1
        pars2[4] = Tsheared2

        return pars1, pars2, pars1_noshear, pars2_noshear


    def _set_shear(self):
        if isinstance(self['shear'], dict):
            shc=self['shear']
            if shc['type']=='normal':
                self._shear_dist=ngmix.priors.RoundGauss2D(shc['mean'][0],
                                                           shc['mean'][1],
                                                           shc['sigma'][0],
                                                           shc['sigma'][1])
            else:
                raise ValueError("only lognormal for now")
            self['shear_is_constant'] = False
        else:
            self._shear_true=self['shear']
            self['shear_is_constant'] = True

    def _get_shear(self):
        """
        if shear is constant, return that, otherwise sample it
        """
        if self['shear_is_constant']:
            shear = self._shear_true
        else:
            shear = self._shear_dist.sample()

        return shear


    def _create_model_set(self):
        """
        crae
        """
        from .sample import DiscreteSampler

        obj_model=self['obj_model']
        if isinstance(obj_model,dict):
            sample_dict = obj_model
        else:
            sample_dict={obj_model: 100}

        self.model_sampler = DiscreteSampler(sample_dict)

        tmod = sample_dict.keys()[0]

    def _set_config(self, sim_conf):
        """
        Check and set the configurations
        """

        sim_conf['do_ring'] = sim_conf.get('do_ring',True)
        if not sim_conf['do_ring']:
            print("    not doing ring")
        self['nsigma_image']=self.get('nsigma_image',NSIGMA_IMAGE)

        self.update(sim_conf)

    def _make_psf(self):
        """
        make the psf gaussian mixture model
        """

        if 'psf_fwhm' in self:
            psf_sigma = self['psf_fwhm']/2.3548200450309493
            self['psf_T'] = 2*psf_sigma**2
            print('psf_T:',self['psf_T'])

        pars=[0.0,
              0.0,
              self['psf_shape'][0],
              self['psf_shape'][1],
              self['psf_T'],
              1.0]

        self.psf_gmix_true=ngmix.gmix.GMixModel(pars, self['psf_model'])
        
    def _get_psf_image(self, dims_pix, cen_arcsec):
        """
        Make the actual image
        """

        self.psf_gmix_true.set_cen(cen_arcsec[0], cen_arcsec[1])
        psf_image=self.psf_gmix_true.make_image(dims_pix,
                                                jacobian=self.jacobian,
                                                nsub=self['nsub'])
        return psf_image

    def _set_pdf(self):
        """
        Set all the priors
        """
        import ngmix
        from ngmix.joint_prior import PriorSimpleSep

        self['pixel_scale'] = self.get('pixel_scale',1.0)

        # we may over-ride below
        self['s2n_for_noise'] = self['s2n']

        # prior separate for all pars
        cen_sigma_arcsec=self['cen_sigma']
        cen_pdf=ngmix.priors.CenPrior(0.0,
                                      0.0,
                                      cen_sigma_arcsec,
                                      cen_sigma_arcsec)

        gtype=self['g_prior_type']
        if gtype=="ba":
            g_pdf_sigma=self['g_prior_sigma']
            g_pdf=ngmix.priors.GPriorBA(g_pdf_sigma)
        elif gtype=="cosmos":
            g_pdf=ngmix.priors.make_gprior_cosmos_sersic(type='erf')
        elif gtype=='great-des':
            g_pdf= ngmix.priors.GPriorGreatDES(pars=self['g_prior_pars'],
                                               gmax=1.0)

        else:
            raise ValueError("only g prior 'ba' for now")

        # T and scatter in linear space
        T            = self['obj_T_mean']
        T_sigma      = self['obj_T_sigma_frac']*T
        counts       = self['obj_counts_mean']
        counts_sigma = self['obj_counts_sigma_frac']*counts


        T_pdf=ngmix.priors.LogNormal(T, T_sigma)
        counts_pdf=ngmix.priors.LogNormal(counts, counts_sigma)

        # for drawing parameters, and after exploration to grab g_pdf and
        # calculate pqr etc.
        self.pdf = PriorSimpleSep(cen_pdf,
                                  g_pdf,
                                  T_pdf,
                                  counts_pdf)

        g_pdf_flat=ngmix.priors.ZDisk2D(1.0)
        self.pdf_gflat = PriorSimpleSep(cen_pdf,
                                        g_pdf_flat,
                                        T_pdf,
                                        counts_pdf)

        self.cen_pdf=cen_pdf
        self.g_pdf=g_pdf
        self.g_pdf_flat=g_pdf_flat
        self.T_pdf=T_pdf
        self.counts_pdf=counts_pdf


    
    def _set_noise(self):
        """
        Find gaussian noise that when added to the image 
        produces the requested s/n.  Use a matched filter.

         sum(pix^2)
        ------------ = S/N^2
          skysig^2

        thus
            
        sum(pix^2)
        ---------- = skysig^2
          (S/N)^2
        """
        
        skysig=self.get('skysig')
        if skysig is not None:
            self['skysig']=skysig
            self['ivar']=1.0/skysig**2
        else:

            imdict=self.get_image_pair_nonoise(random=False)
            im=imdict['im1']['obs0'].image
            skysig2 = (im**2).sum()/self['s2n_for_noise']**2
            skysig = numpy.sqrt(skysig2)


            self['skysig']=skysig
            self['ivar']=1.0/skysig**2

            imn, noise_im =self._add_noise(im)
            s2n_numer = (imn*im*self['ivar']).sum()
            s2n_denom = numpy.sqrt( (im**2*self['ivar']).sum() )
            s2n_check = s2n_numer/s2n_denom


    def _get_model_s2n(self, im):
        s2n = numpy.sqrt( (im**2).sum() )/self['skysig']
        return s2n


    def _add_noise(self, im):
        """
        Add gaussian random noise
        """
        nim = self['skysig']*randn(im.size).reshape(im.shape)
        return im + nim, nim


    def _get_dims_cen(self, T):
        """
        Based on T, get the required dimensions and a center

        Should send T+Tpsf
        """

        sigma_pix=numpy.sqrt(T/2.)/self['pixel_scale']
        dims_pix = array( [2.*sigma_pix*self['nsigma_image']]*2 )

        cen_pix = array( [(dims_pix[0]-1.)/2.]*2 )

        return dims_pix, cen_pix


def srandu(num=None):
    """
    Generate random numbers in the symmetric distribution [-1,1]
    """
    return 2*(numpy.random.random(num)-0.5)



def seed_global_devrand():
    """
    Seed the "global" random number generator
    """
    seed = get_devrand_uint()
    numpy.random.seed(seed)

def get_random_state_devrand(seed=None):
    """
    get the random state.  If seed is not sent, 
    get the eed the numpy random state from /dev/random
    """

    if seed is None:
        seed = get_devrand_uint()
    return numpy.random.mtrand.RandomState(seed=seed)

def get_devrand_uint():
    """
    Read 4 bytes from /dev/random and convert to an
    unsigned int. The returned value is a normal
    python int
    """
    import struct
    four_bytes = get_random_bytes(4)
    # I is unsigned int
    tup = struct.unpack("I", four_bytes)
    val = tup[0]
    return val


def get_random_bytes(nbytes):
    """
    Get the specified number of bytes from /dev/random
    """
    fd = os.open("/dev/random",os.O_RDONLY)
    thebytes=os.read(fd, 4)
    os.close(fd)

    return thebytes
