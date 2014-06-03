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
from numpy import array, zeros, log, log10, exp, sqrt
from numpy.random import random as randu
from numpy.random import randn

import ngmix

# region over which to render images and calculate likelihoods
NSIGMA_RENDER=5.0

# minutes
DEFAULT_CHECKPOINTS=[30,60,90,110]
#DEFAULT_CHECKPOINTS=[30,90]
#DEFAULT_CHECKPOINTS=[0,100]

class TryAgainError(Exception):
    def __init__(self, message):

        # Call the base class constructor with the parameters it needs
        Exception.__init__(self, message)


class NGMixSim(dict):
    def __init__(self, sim_conf, run_conf, s2n, npairs, **keys):
        """
        Simulate and fit the requested number of pairs at
        the specified s/n
        """

        # seed a new MT random generator from devrand
        # and return the object
        self.random_state=get_random_state_devrand()
        # also seed the global random number generator from devrand
        seed_global_devrand()

        self.set_config(sim_conf, run_conf)
        self.update(self.conf)
        self.update(keys)

        self.shear=self.simc['shear']
        self.nsub=self.simc['nsub']
        self.nsigma_render=self.simc.get('nsigma_render',NSIGMA_RENDER)

        self.check_pqr_shear()

        # this might be mean for correspond to a flux threshold, 
        # depending on the type of run
        self.s2n=s2n
        self.npairs=npairs

        self.obj_model=self.simc['obj_model']
        self.fit_model=self['fit_model']

        self.true_npars = ngmix.gmix.get_model_npars(self.obj_model)

        if self.fit_model is not None:
            self.npars=ngmix.gmix.get_model_npars(self.fit_model)
            print("npars:",self.npars)

        self.make_plots=keys.get('make_plots',False)
        self.separate=False
        self.plot_base=keys.get('plot_base',None)

        self.setup_checkpoints(**keys)

        self.verbose=run_conf.get('verbose',True)

        if self.data is None:
            self.make_struct()

        self.set_prior()
        self.make_psf()

        self.set_noise()

        if self.verbose:
            pprint.pprint(self)
            pprint.pprint(self.simc)

    def run_sim(self):
        """
        Run the simulation, fitting psf and all pairs
        """

        self.start_timer()

        i=0
        npairs=self.npairs
        for ipair in xrange(npairs):
            if self.verbose:
                print('%s/%s' % (ipair+1,npairs) )

            self.ipair=ipair
            if self.data['processed'][i]:
                i += 2 # skip the pair
            else:
                while True:
                    try:
                        reslist=self.process_pair()
                        break
                    except TryAgainError as err:
                        if self.verbose:
                            print(str(err))

                self.copy_to_output(reslist[0], i)
                i += 1
                self.copy_to_output(reslist[1], i)
                i += 1

            self.set_elapsed_time()
            self.try_checkpoint()

        self.set_elapsed_time()
        if self.verbose:
            print('time minutes:',self.tm_minutes)
            print('time per pair sec:',self.tm/npairs)
            print('time per image sec:',self.tm/(2*npairs))



    def set_config(self, sim_conf, run_conf):
        """
        Check and set the configurations
        """
        if sim_conf['name'] != run_conf['sim']:
            err="sim name in run config '%s' doesn't match sim name '%s'"
            raise ValueError(err % (run_conf['sim'],sim_conf['name']))

        self.simc=sim_conf
        self.conf=run_conf

    def check_pqr_shear(self):
        """
        If we are doing pqr need to check what shear we expand
        about
        """

        expand_shear_true=self.get('expand_shear_true',None)
        if expand_shear_true:
            self.shear_expand=array(self.shear)
            print('nsim: expanding about shear:',self.shear_expand)
        else:
            self.shear_expand=array([0.0,0.0])

    def get_data(self):
        """
        Get a ref to the data array with the fit results
        """
        return self.data


    def start_timer(self):
        """
        Set the elapsed time so far
        """
        self.tm0 = time.time()

    def set_elapsed_time(self):
        """
        Set the elapsed time so far
        """

        self.tm = time.time()-self.tm0
        self.tm_minutes = self.tm/60.0

    def process_pair(self):
        """
        Create a simulated image pair and perform the fit
        """

        imdicts = self.get_noisy_image_pair()

        print(imdicts['im1']['obs'].image.shape)

        reslist=[]
        for key in ['im1','im2']:

            imd = imdicts[key]

            obs = imd['obs']

            # note obs.psf is another Observation.  the gmix for
            # psf will be set
            self.fit_psf(obs.psf)

            fitter=self.fit_galaxy(imd)
            res=fitter.get_result()
            if res['flags'] != 0:
                raise TryAgainError("failed at %s flags %s" % (key,res['flags']))

            res['pars_true'] = imd['pars']
            res['s2n_true'] = imd['s2n']
            self.print_res(res)

            reslist.append(res)

            if self.make_plots:
                self.do_make_plots(fitter,key)

        return reslist

    def do_make_plots(self, fitter, key):
        """
        Write a plot file of the trials
        """
        import biggles
        biggles.configure('default','fontsize_min',1.0)

        width,height=1100,1100

        pdict=fitter.make_plots(title=self.fit_model,
                                separate=self.separate,
                                weights=self._weights)

        if not self.separate:

            trials_pname=self.plot_base+'-%06d-%s-trials.png' % (self.ipair,key)
            print(trials_pname)
            p=pdict['trials']
            p.write_img(width,height,trials_pname)

            if 'wtrials' in pdict:
                wp=pdict['wtrials']
                wtrials_pname=\
                    self.plot_base+'-%06d-%s-wtrials.png' % (self.ipair,key)
                print(wtrials_pname)
                wp.write_img(width,height,wtrials_pname)
        else:
            burn_plt, hist_plt=pdict['trials']
            burn_pname=self.plot_base+'-%06d-%s-burn.png' % (self.ipair,key)
            hist_pname=self.plot_base+'-%06d-%s-hist.png' % (self.ipair,key)

            print(burn_pname)
            burn_plt.write_img(width,height,burn_pname)
            print(hist_pname)
            hist_plt.write_img(width,height,hist_pname)

            if 'wtrials' in pdict:
                wburn_plt, whist_plt=pdict['wtrials']
                whist_pname=self.plot_base+'-%06d-%s-whist.png' % (self.ipair,key)

                print(whist_pname)
                whist_plt.write_img(width,height,whist_pname)


    def fit_galaxy(self, imdict):
        """
        Fit according to the requested method
        """

        fitter_type=self['fitter']
        #if fitter_type == 'mcmc':
        if 'mcmc' in fitter_type:
            fitter = self.fit_galaxy_mcmc(imdict)
        elif 'mh' in fitter_type:
            fitter = self.fit_galaxy_mh(imdict)
        elif fitter_type == 'lm':
            fitter = self.fit_galaxy_lm(imdict,
                                        prior_type=self['prior_type_during'],
                                        ntry=self['lm_ntry'])

        else:
            raise ValueError("bad fitter type: '%s'" % fitter_type)

        return fitter

    def fit_galaxy_mcmc(self, imdict):
        """
        Fit the model to the galaxy
        """
        fitter=self.run_simple_mcmc_fitter(imdict)
        self._add_mcmc_stats(fitter)
        return fitter

    def fit_galaxy_mh(self, imdict):
        """
        Fit the model to the galaxy
        """
        fitter=self.run_simple_mh_fitter(imdict)
        self._add_mcmc_stats(fitter)
        return fitter


    def _add_mcmc_stats(self, fitter):
        """
        Calculate some stats

        The result dict internal to the fitter is modified to include
        g_sens and P,Q,R

        could we add this as a stand-alone function to ngmix.fitting?
        """

        trials = fitter.get_trials()
        g=trials[:,2:2+2]

        # this is the full prior
        g_prior=self.prior.g_prior
        weights = g_prior.get_prob_array2d(g[:,0], g[:,1])

        # keep for later if we want to make plots
        self._weights=weights

        fitter.calc_result(weights=weights)

        # we are going to mutate the result dict owned by the fitter
        res=fitter.get_result()

        ls=ngmix.lensfit.LensfitSensitivity(g, g_prior)
        g_sens = ls.get_g_sens()
        g_mean = ls.get_g_mean()
        nuse = ls.get_nuse()

        pqrobj=ngmix.pqr.PQR(g, g_prior,
                             shear_expand=self.shear_expand)
        P,Q,R = pqrobj.get_pqr()

        # this nuse should be the same for both lensfit and pqr
        res['nuse'] = nuse
        res['g_sens'] = g_sens
        res['P']=P
        res['Q']=Q
        res['R']=R


    def get_guess(self, imdict, n=1):
        """
        Get a full guess, nwalkers x npars
        """

        guess_type=self['guess_type']

        if guess_type=='draw_truth':
            print('    * guessing randomized truth')
            full_guess=self.get_guess_from_pars(imdict['pars'], n=n)
        elif guess_type=='draw_priors':
            full_guess=self.get_guess_draw_priors(n=n)
        elif guess_type=='draw_maxlike':
            full_guess,perr=self.get_guess_draw_maxlike(imdict, n=n)
        else:
            raise ValueError("bad guess type: '%s'" % guess_type)
        
        return full_guess


    def run_simple_mcmc_fitter(self, imdict):
        """
        simple gauss,exp,dev
        """
        from ngmix.fitting import MCMCSimple

        guess=self.get_guess(imdict, n=self['nwalkers'])

        obs=imdict['obs']

        fitter=MCMCSimple(obs,
                          self.fit_model,

                          prior=self.prior_gflat, # no prior during

                          nwalkers=self['nwalkers'],
                          mca_a=self['mca_a'],
                          random_state=self.random_state)

        pos=fitter.run_mcmc(guess,self['burnin'])
        pos=fitter.run_mcmc(pos,self['nstep'])

        return fitter

    def run_simple_mh_fitter(self, imdict):
        """
        simple gauss,exp,dev
        """
        from ngmix.fitting import MHSimple

        mess="for mh guess should be maxlike"
        assert (self['guess_type']=="draw_maxlike"),mess

        guess,perr=self.get_guess_draw_maxlike(imdict, n=1)

        step_sizes = 0.5*perr

        obs=imdict['obs']

        fitter=MHSimple(obs,
                        self.fit_model,
                        step_sizes,

                        prior=self.prior_gflat)

        pos=fitter.run_mcmc(guess,self['burnin'])
        pos=fitter.run_mcmc(pos,self['nstep'])

        return fitter


    def get_guess_draw_priors(self, n=1):
        """
        Get a guess drawn from the priors

        assume simple for now
        """
        print("drawing guess from priors")
        guess=self.prior.sample(n)

        if n==1:
            guess=guess[0,:]

        return guess

    def get_guess_draw_maxlike(self, imdict, n=1):
        """
        Get the maximum likelihood fit and draw from that using
        width from the fit

        start lm from priors to make it a challenge.  always use full prior to
        make sure we get a good estimate of the covariance matrix (without
        priors the covar can blow up)

        """

        print("drawing guess from maxlike")
        ntry=self['lm_ntry']
        for i in xrange(ntry):

            #lm_guess=self.get_guess_from_pars(imdict['pars'], n=1)
            lm_guess=self.get_guess_draw_priors(n=1)
            fitter=self.fit_galaxy_lm(imdict,
                                      guess=lm_guess,
                                      prior_type='full',
                                      ntry=1)
            res=fitter.get_result()
            if res['flags']==0:
                break

        if res['flags'] != 0:
            raise TryAgainError("failed fit LM after %s tries" % ntry)

        pars=res['pars']
        perr=res['pars_err']
        ngmix.fitting.print_pars(pars, front='        lmpars: ')
        ngmix.fitting.print_pars(perr, front='        lmperr: ')

        # now scatter it around a bit
        guess=self.get_guess_from_pars(pars, n=n, width=perr)

        return guess, perr

    def get_guess_from_pars(self, pars, n=1, width=None):
        """
        Get a guess centered on the input pars

        pars are in log everywhere
        """

        npars=pars.size

        if width is None:
            width = pars*0 + 0.01
        else:
            if len(width) != len(pars):
                raise ValueError("width not same size as pars")

        guess=numpy.zeros( (n, pars.size) )

        guess[:,0] = width[0]*srandu(n)
        guess[:,1] = width[1]*srandu(n)

        guess_shape=self.get_shape_guess(pars[2],pars[3],n,width=width[2:2+2])
        guess[:,2]=guess_shape[:,0]
        guess[:,3]=guess_shape[:,1]

        # we add to log pars!
        for i in xrange(4,npars):
            guess[:,i] = pars[i] + width[i]*srandu(n)

        if n==1:
            guess=guess[0,:]

        return guess

    def get_shape_guess(self, g1, g2, n, width=[0.01,0.01]):
        """
        Get guess, making sure in range
        """
        from ngmix.gexceptions import GMixRangeError

        guess=numpy.zeros( (n, 2) )
        shape=ngmix.Shape(g1, g2)

        for i in xrange(n):

            while True:
                try:
                    g1_offset = width[0]*srandu()
                    g2_offset = width[1]*srandu()
                    shape_new=shape.copy()
                    shape_new.shear(g1_offset, g2_offset)
                    break
                except GMixRangeError:
                    pass

            guess[i,0] = shape_new.g1
            guess[i,1] = shape_new.g2

        return guess


    def fit_galaxy_lm(self, imdict, prior_type='full', guess=None, ntry=1):
        """
        Fit the model to the galaxy

        we send some keywords so behavior can be different if this is a guess
        getter

        """
        from ngmix.fitting import LMSimple

        obs=imdict['obs']

        # no prior for now
        if prior_type=='full':
            print("using full prior for lm")
            prior=self.prior
        elif prior_type=='gflat':
            print("using flat g prior for lm")
            prior=self.prior_gflat
        elif prior_type==None:
            prior=None
        else:
            raise ValueError("bad prior type during: %s" % ptd)

        for i in xrange(ntry):
            if guess is None or i > 0:
                guess=self.get_lm_guess(imdict)

            fitter=LMSimple(obs,
                            self.fit_model,

                            prior=prior,

                            lm_pars=self['lm_pars'])

            fitter.run_lm(guess)
            res=fitter.get_result()
            if res['flags']==0:
                break

        res['ntry']=i+1
        return fitter


    def get_lm_guess(self, imdict):
        """
        Get a guess for the LM fitter
        """
        guess_type=self['guess_type']
        if guess_type=='truth':
            guess=imdict['pars']
        elif guess_type=='truth_random':
            guess=self.get_guess_from_pars(imdict['pars'],n=1)
        elif guess_type=='draw_priors':
            guess=self.get_guess_draw_priors(n=1)
        else:
            raise ValueError("bad guess type: '%s'" % guess_type)
        return guess

    def print_res(self,res):
        """
        print some stats
        """

        if 'arate' in res:
            print('    arate:',res['arate'],'s2n_w:',
                    res['s2n_w'],'nuse:',res['nuse'])
        elif 'nfev' in res:
            print('    nfev:',res['nfev'],'s2n_w:',res['s2n_w'])

        ngmix.fitting.print_pars(res['pars_true'], front='    true: ')
        ngmix.fitting.print_pars(res['pars'],front='    pars: ')
        ngmix.fitting.print_pars(res['pars_err'],front='    perr: ')

    def make_psf(self):
        """
        make the psf gaussian mixture model
        """

        if 'psf_fwhm' in self.simc:
            psf_sigma = self.simc['psf_fwhm']/2.3548200450309493
            self.simc['psf_T'] = 2*psf_sigma**2
            print('psf_T:',self.simc['psf_T'])

        pars=[0.0,
              0.0,
              self.simc['psf_shape'][0],
              self.simc['psf_shape'][1],
              self.simc['psf_T'],
              1.0]

        # note not log pars
        self.psf_gmix_true=ngmix.gmix.GMixModel(pars, self.simc['psf_model'])
        
    def get_psf_image(self, dims_pix, cen_arcsec):
        """
        Make the actual image
        """

        self.psf_gmix_true.set_cen(cen_arcsec[0], cen_arcsec[1])
        psf_image=self.psf_gmix_true.make_image(dims_pix,
                                                jacobian=self.jacobian,
                                                nsub=self.nsub)
        return psf_image

    def fit_psf(self, obs):
        """
        Fit the pixelized psf to a model

        will set_gmix
        """
        from ngmix.observation import Observation
        from ngmix.gexceptions import GMixMaxIterEM

        print('    fitting psf')

        image=obs.image
        imsky,sky=ngmix.em.prep_image(image)

        sky_obs = Observation(imsky, jacobian=obs.jacobian)

        em=ngmix.em.GMixEM(sky_obs)

        tol=self.get('psf_tol',1.0e-6)
        maxiter=self.get('psf_maxiter',30000)

        while True:
            guess=self.get_psf_guess()
            print('    psf guess:')
            print(guess)
            try:
                em.go(guess, sky, tol=tol,maxiter=maxiter)
                break
            except GMixMaxIterEM as e:
                print(str(e))
                print('re-trying')

        psf_gmix_fit=em.get_gmix()

        #print('psf fit:')
        #print(psf_gmix_fit)

        obs.set_gmix(psf_gmix_fit)

    def get_psf_guess(self):
        """
        Get the starting guess.

        If psf_ngauss is not sent, use the "truth" as the guess.  This might
        not be as good a fit if for example an extra gaussian is needed to deal
        with pixelization.  Probably more important for galaxies ~size of the
        """
        import ngmix

        ngauss=self.get('psf_ngauss',None)

        pt=self.psf_gmix_true
        ngauss_true=len(pt)
        if ngauss is None or ngauss==ngauss_true:
            # we can just use the "truth" as a guess
            guess=pt.copy()
        else:
            # more generic guess
            T = pt.get_T()
            cen = pt.get_cen()

            pars=numpy.zeros(ngauss*6)

            if ngauss==2:
                p1 = 0.6*(1.0 + 0.05*srandu() )
                p2 = 0.4*(1.0 + 0.05*srandu() )
                T1 = T*0.4*(1.0 + 0.1*srandu() )
                T2 = T*0.6*(1.0 + 0.1*srandu() )

                pars[0] = p1
                pars[1] = cen[0] + 0.1*srandu()
                pars[2] = cen[1] + 0.1*srandu()
                pars[3] = T1/2.
                pars[4] = 0.0
                pars[5] = T1/2.

                pars[6] = p2
                pars[7] = cen[0] + 0.1*srandu()
                pars[8] = cen[1] + 0.1*srandu()
                pars[9] = T2/2.
                pars[10] = 0.0
                pars[11] = T2/2.

            else: 
                raise ValueError("support ngauss > 2")
            
            # not log pars
            guess=ngmix.gmix.GMix(pars=pars)

        return guess

    def set_prior(self):
        """
        Set all the priors
        """
        import ngmix

        self.pixel_scale = self.simc.get('pixel_scale',1.0)

        # we may over-ride below
        self.s2n_for_noise = self.s2n

        simc=self.simc

        if simc['prior_type'] == "separate":
            from ngmix.joint_prior import PriorSimpleSep
            # prior separate for all pars
            cen_sigma_arcsec=simc['cen_sigma']
            cen_prior=ngmix.priors.CenPrior(0.0,
                                            0.0,
                                            cen_sigma_arcsec,
                                            cen_sigma_arcsec)

            gtype=simc['g_prior_type']
            if gtype=="ba":
                g_prior_sigma=simc['g_prior_sigma']
                g_prior=ngmix.priors.GPriorBA(g_prior_sigma)
            elif gtype=="cosmos":
                g_prior=ngmix.priors.make_gprior_cosmos_sersic(type='erf')
            else:
                raise ValueError("only g prior 'ba' for now")

            g_prior_flat=ngmix.priors.ZDisk2D(1.0)

            # T and scatter in linear space, convert to log
            T            = simc['obj_T_mean']
            T_sigma      = simc['obj_T_sigma_frac']*T
            counts       = simc['obj_counts_mean']
            counts_sigma = simc['obj_counts_sigma_frac']*counts

            logT_mean, logT_sigma=ngmix.priors.lognorm_convert(T,T_sigma,base=10.0)
            logcounts_mean, logcounts_sigma=ngmix.priors.lognorm_convert(counts,counts_sigma,base=10.0)

            T_prior=ngmix.priors.Normal(logT_mean, logT_sigma)
            counts_prior=ngmix.priors.Normal(logcounts_mean, logcounts_sigma)

            # for drawing parameters, and after exploration to grab g_prior and calculate
            # pqr etc.
            self.prior = PriorSimpleSep(cen_prior,
                                        g_prior,
                                        T_prior,
                                        counts_prior)

            # for the exploration, for which we do not apply g prior during
            self.prior_gflat = PriorSimpleSep(cen_prior,
                                              g_prior_flat,
                                              T_prior,
                                              counts_prior)

        else:
            raise ValueError("bad prior type: %s" % simc['prior_type'])

    
    def set_noise(self):
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
        
        skysig=self.simc.get('skysig')
        if skysig is not None:
            self.skysig=skysig
            self.ivar=1.0/skysig**2
        else:

            imdict=self.get_image_pair(random=False)
            im=imdict['im1']['obs0'].image
            skysig2 = (im**2).sum()/self.s2n_for_noise**2
            skysig = numpy.sqrt(skysig2)


            self.skysig=skysig
            self.ivar=1.0/skysig**2

            imn=self.add_noise(im)
            s2n_numer = (imn*im*self.ivar).sum()
            s2n_denom = numpy.sqrt( (im**2*self.ivar).sum() )
            s2n_check = s2n_numer/s2n_denom


    def get_model_s2n(self, im):
        s2n = numpy.sqrt( (im**2).sum() )/self.skysig
        return s2n

    def get_noisy_image_pair(self):
        """
        Get an image pair, with noise added
        """
        from ngmix.observation import Observation

        imdict=self.get_image_pair()

        obs1_0 = imdict['im1']['obs0']
        obs2_0 = imdict['im2']['obs0']

        im1_0 = obs1_0.image
        im2_0 = obs2_0.image
        im1=self.add_noise(im1_0)
        im2=self.add_noise(im2_0)

        wt=numpy.zeros(im1_0.shape) + self.ivar

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

        im1_s2n = self.get_model_s2n(im1_0)
        im2_s2n = self.get_model_s2n(im2_0)

        imdict['im1']['s2n'] = im1_s2n
        imdict['im2']['s2n'] = im2_s2n

        return imdict

    def add_noise(self, im):
        """
        Add gaussian random noise
        """

        return im + self.skysig*randn(im.size).reshape(im.shape)

    def get_image_pair(self, random=True):
        """
        get a model image

        If random is True, use draw random values from the priors.
        Otherwise use the mean of the priors

        """
        import ngmix
        from ngmix.observation import Observation

        # pars are in log space
        pars1, pars2 = self.get_pair_pars(random=random)

        gm1_pre=ngmix.gmix.GMixModel(pars1, self.obj_model, logpars=True)
        gm2_pre=ngmix.gmix.GMixModel(pars2, self.obj_model, logpars=True)

        gm1  = gm1_pre.convolve(self.psf_gmix_true)
        gm2  = gm2_pre.convolve(self.psf_gmix_true)

        T = gm1.get_T()
        dims_pix, cen0_pix = self.get_dims_cen(T)

        # conversion between pixels and sky in arcsec
        self.jacobian=ngmix.jacobian.Jacobian(cen0_pix[0],
                                              cen0_pix[1],
                                              self.pixel_scale,
                                              0.0,
                                              0.0,
                                              self.pixel_scale)
        psf_cen=pars1[0:0+2]
        psf_image=self.get_psf_image(dims_pix, psf_cen)

        psf_obs = Observation(psf_image, jacobian=self.jacobian)

        nsub = self.nsub
        im1=gm1.make_image(dims_pix, nsub=nsub, jacobian=self.jacobian)
        im2=gm2.make_image(dims_pix, nsub=nsub, jacobian=self.jacobian)

        obs1 = Observation(im1, psf=psf_obs, jacobian=self.jacobian)
        obs2 = Observation(im2, psf=psf_obs, jacobian=self.jacobian)
        
        out={'im1':{'pars':pars1,'gm_pre':gm1_pre,'gm':gm1,'obs0':obs1},
             'im2':{'pars':pars2,'gm_pre':gm2_pre,'gm':gm2,'obs0':obs2}}
        return out

    def get_pair_pars(self, random=True):
        """
        Get pair parameters

        if not random, then the mean pars are used, except for cen and g1,g2
        which are zero

        Note the pars are in log space for T,F
        """
        import ngmix
        from numpy import pi

        if random:
            pars1 = self.prior.sample()

            # use same everything but rotated 90 degrees
            pars2=pars1.copy()

            g1=pars1[2]
            g2=pars1[3]
        
            shape2=ngmix.shape.Shape(g1,g2)
            shape2.rotate(pi/2.0)

            pars2[2]=shape2.g1
            pars2[3]=shape2.g2

        else:
            samples=self.prior.sample(10000)
            pars1 = samples.mean(axis=0)

            pars1[0:0+4] = 0.0
            pars2=pars1.copy()

        shape1=ngmix.shape.Shape(pars1[2],pars1[3])
        shape2=ngmix.shape.Shape(pars2[2],pars2[3])

        shear=self.shear
        shape1.shear(shear[0], shear[1])
        shape2.shear(shear[0], shear[1])

        pars1[2]=shape1.g1
        pars1[3]=shape1.g2
        pars2[2]=shape2.g1
        pars2[3]=shape2.g2

        return pars1, pars2

    def get_dims_cen(self, T):
        """
        Based on T, get the required dimensions and a center
        """

        sigma_pix=numpy.sqrt(T/2.)/self.pixel_scale

        dims_pix = array( [2.*sigma_pix*self.nsigma_render]*2 )
        cen_pix = array( [(dims_pix[0]-1.)/2.]*2 )

        return dims_pix, cen_pix

    def setup_checkpoints(self, **keys):
        """
        Set up checkpoint times, file, and sent data
        """

        self.checkpoints     = keys.get('checkpoints',DEFAULT_CHECKPOINTS)
        self.n_checkpoint    = len(self.checkpoints)
        self.checkpointed    = [False]*self.n_checkpoint

        self.checkpoint_file=keys.get('checkpoint_file',None)
        self.set_checkpoint_data(**keys)

        if self.checkpoint_file is not None:
            self.do_checkpoint=True
        else:
            self.do_checkpoint=False


    def set_checkpoint_data(self, **keys):
        """
        Look for checkpoint data, file etc.
        """
        self.data=None

        checkpoint_data=keys.get('checkpoint_data',None)
        if checkpoint_data is not None:
            self.data=checkpoint_data


    def try_checkpoint(self):
        """
        If we should make a checkpoint, do so
        """

        should_checkpoint, icheck = self.should_checkpoint()

        if should_checkpoint:
            self.write_checkpoint()
            self.checkpointed[icheck]=True


    def should_checkpoint(self):
        """
        Should we write a checkpoint file?
        """

        should_checkpoint=False
        icheck=-1

        if self.do_checkpoint:
            for i in xrange(self.n_checkpoint):

                checkpoint=self.checkpoints[i]
                checkpointed=self.checkpointed[i]

                if self.tm_minutes > checkpoint and not checkpointed:
                    should_checkpoint=True
                    icheck=i

        return should_checkpoint, icheck


    def write_checkpoint(self):
        """
        Write the checkpoint file

        The file is written to the cwd and if this is not the final
        destination, an attempt is made to move it there.  This may fail and if
        so a message is printed.

        """

        print('checkpointing at',self.tm_minutes,'minutes')
        success=write_fits(self.checkpoint_file, self.data)


    def copy_to_output(self, res, i):
        """
        Copy results into the output
        """
        d=self.data
        d['processed'][i] = 1

        d['s2n_true'][i] = res['s2n_true']
        d['pars_true'][i,:] = res['pars_true']

        d['pars'][i,:] = res['pars']
        d['pcov'][i,:,:] = res['pars_cov']

        d['g'][i,:] = res['g']
        d['g_cov'][i,:,:] = res['g_cov']

        d['s2n_w'][i] = res['s2n_w']

        if 'P' in res:
            d['arate'][i] = res['arate']
            d['g_sens'][i,:] = res['g_sens']

            d['P'][i] = res['P']
            d['Q'][i,:] = res['Q']
            d['R'][i,:,:] = res['R']

            d['nuse'][i] = res['nuse']
        else:
            d['nfev'][i] = res['nfev']
            # set outside of fitter
            d['ntry'][i] = res['ntry']

    def make_struct(self):
        """
        Make the output array
        """
        npars=self.npars

        dt=[('processed','i2'),
            ('s2n_true','f8'),
            ('pars_true','f8',self.true_npars),
            ('pars','f8',npars),
            ('pcov','f8',(npars,npars)),
            ('s2n_w','f8'),
            ('arate','f8')]

        if 'lm' in self['fitter']:
            dt += [('nfev','i4'),
                   ('ntry','i4'),
                   ('g','f8',2),
                   ('g_cov','f8',(2,2))]
        else:
            dt += [('g','f8',2),
                   ('g_cov','f8',(2,2)),
                   ('g_sens','f8',2),

                   ('P','f8'),
                   ('Q','f8',2),
                   ('R','f8',(2,2)),
                   ('nuse','i4')]
        self.data=numpy.zeros(self.npairs*2, dtype=dt)




class NGMixSimJointSimpleLinPars(NGMixSim):
    """
    Simple model with joint prior on g1,g2,T,Flux
    """
    def __init__(self, sim_conf, run_conf, npairs, **keys):
        s2n=-1
        super(NGMixSimJointSimpleLinPars,self).__init__(sim_conf, run_conf, s2n, npairs, **keys)

    def fit_galaxy(self, imdict):
        """
        Fit simple model with joint prior
        """
        from ngmix.fitting import MCMCSimpleJointLinPars

        nwalkers = self['nwalkers']

        full_guess=self.get_guess(imdict, n=nwalkers)

        fitter=MCMCSimpleJointLinPars(imdict['image'],
                                      imdict['wt'],
                                      self.jacobian,
                                      self['fit_model'],

                                      cen_prior=self.cen_prior,
                                      joint_prior=self.joint_prior,

                                      full_guess=full_guess,

                                      shear_expand=self.shear_expand,

                                      psf=self.psf_gmix_fit,
                                      nwalkers=nwalkers,
                                      nstep=self['nstep'],
                                      burnin=self['burnin'],
                                      mca_a=self['mca_a'],

                                      prior_during=self['prior_during'],

                                      random_state=self.random_state,
                                      do_pqr=True)
        fitter.go()
        res=fitter.get_result()

        return fitter


    def get_pair_pars(self, **keys):
        """
        Get pair parameters
        """
        from numpy import pi
        import ngmix

        # centers zero
        pars1=numpy.zeros(6)

        # [g,T,Flux]
        sample=self.joint_prior.sample2d()
        pars1[4] = sample[2]
        pars1[5] = sample[3]

        g1,g2=sample[0],sample[1]
        shape1=ngmix.shape.Shape(g1,g2)

        shape2=shape1.copy()
        shape2.rotate(pi/2.0)

        shear=self.shear
        shape1.shear(shear[0], shear[1])
        shape2.shear(shear[0], shear[1])

        pars1[2] = shape1.g1
        pars1[3] = shape1.g2

        pars2=pars1.copy()
        pars2[2] = shape2.g1
        pars2[3] = shape2.g2

        cen_offset_arcsec=array( self.cen_prior.sample() )
        return pars1, pars2, cen_offset_arcsec

    def set_prior(self):
        """
        Set all the priors
        """
        import great3
        import ngmix


        self.pixel_scale = self.simc.get('pixel_scale',1.0)

        joint_prior_type = self.simc['joint_prior_type']
        self.joint_prior = \
            great3.joint_prior.make_joint_prior_simple(type=joint_prior_type)

        cen_sigma_arcsec=self.simc['cen_sigma']
        self.cen_prior=ngmix.priors.CenPrior(0.0,
                                             0.0,
                                             cen_sigma_arcsec,
                                             cen_sigma_arcsec)




class NGMixSimJointSimpleHybrid(NGMixSim):
    """
    Simple model with joint prior on g1,g2,T,Flux
    """
    def __init__(self, sim_conf, run_conf, npairs, **keys):
        s2n=-1
        super(NGMixSimJointSimpleHybrid,self).__init__(sim_conf, run_conf, s2n, npairs, **keys)

    def fit_galaxy(self, imdict):
        """
        Fit simple model with joint prior
        """
        from ngmix.fitting import MCMCSimpleJointHybrid

        nwalkers = self['nwalkers']

        full_guess=self.get_guess(imdict['pars'], n=nwalkers)

        fitter=MCMCSimpleJointHybrid(imdict['image'],
                                     imdict['wt'],
                                     self.jacobian,
                                     self['fit_model'],

                                     cen_prior=self.cen_prior,
                                     joint_prior=self.joint_prior,

                                     full_guess=full_guess,

                                     shear_expand=self.shear_expand,

                                     psf=self.psf_gmix_fit,
                                     nwalkers=nwalkers,
                                     nstep=self['nstep'],
                                     burnin=self['burnin'],
                                     mca_a=self['mca_a'],

                                     prior_during=self['prior_during'],

                                     random_state=self.random_state,
                                     do_pqr=True)
        fitter.go()
        res=fitter.get_result()

        return fitter


    def get_pair_pars(self, **keys):
        """
        Get pair parameters
        """
        from numpy import pi
        import ngmix

        # centers zero
        pars1=numpy.zeros(6)

        g1a,g2a = self.joint_prior.g_prior.sample2d(1)

        sample=self.joint_prior.sample()
        logT=sample[0]
        logF=sample[1]

        pars1[4] = 10.0**logT
        pars1[5] = 10.0**logF

        shape1=ngmix.shape.Shape(g1a[0],g2a[0])

        shape2=shape1.copy()
        shape2.rotate(pi/2.0)

        shear=self.shear
        shape1.shear(shear[0], shear[1])
        shape2.shear(shear[0], shear[1])

        pars1[2] = shape1.g1
        pars1[3] = shape1.g2

        pars2=pars1.copy()
        pars2[2] = shape2.g1
        pars2[3] = shape2.g2

        print("log pars:",logT,logF)

        cen_offset_arcsec=array( self.cen_prior.sample() )
        return pars1, pars2, cen_offset_arcsec

    def set_prior(self):
        """
        Set all the priors
        """
        import great3
        import ngmix


        self.pixel_scale = self.simc.get('pixel_scale',1.0)

        joint_prior_type = self.simc['joint_prior_type']
        self.joint_prior = \
            great3.joint_prior.make_joint_prior_simple(type=joint_prior_type)

        cen_sigma_arcsec=self.simc['cen_sigma']
        self.cen_prior=ngmix.priors.CenPrior(0.0,
                                             0.0,
                                             cen_sigma_arcsec,
                                             cen_sigma_arcsec)


    def get_guess(self, pars, n=1, width=None):
        """
        Get a guess centered on the truth

        width is relative for T and counts
        """

        if width is None:
            width = pars*0 + 0.01
        else:
            if len(width) != len(pars):
                raise ValueError("width not same size as pars")

        guess=numpy.zeros( (n, pars.size) )

        guess[:,0] = width[0]*srandu(n)
        guess[:,1] = width[1]*srandu(n)
        guess_shape=self.get_shape_guess(pars[2],pars[3],n,width=width[2:2+2])

        guess[:,2]=guess_shape[:,0]
        guess[:,3]=guess_shape[:,1]

        logT=log10(pars[4])
        logF=log10(pars[5])

        guess[:,4] = logT*(1.0 + width[4]*srandu(n))
        guess[:,5] = logF*(1.0 + width[5]*srandu(n))


        return guess



class NGMixSimJointSimpleLogPars(NGMixSim):
    """
    Simple model with joint prior on g1,g2,T,Flux
    """
    def __init__(self, sim_conf, run_conf, npairs, **keys):
        s2n=-1
        super(NGMixSimJointSimple,self).__init__(sim_conf, run_conf, s2n, npairs, **keys)

    def fit_galaxy(self, imdict):
        """
        Fit simple model with joint prior
        """
        from ngmix.fitting import MCMCSimpleJoint

        nwalkers = self['nwalkers']

        full_guess=self.get_guess(imdict['pars'], n=nwalkers)

        fitter=MCMCSimpleJoint(imdict['image'],
                               imdict['wt'],
                               self.jacobian,
                               self['fit_model'],

                               cen_prior=self.cen_prior,
                               joint_prior=self.joint_prior,

                               full_guess=full_guess,

                               shear_expand=self.shear_expand,

                               psf=self.psf_gmix_fit,
                               nwalkers=nwalkers,
                               nstep=self['nstep'],
                               burnin=self['burnin'],
                               mca_a=self['mca_a'],

                               prior_during=self['prior_during'],

                               random_state=self.random_state,
                               do_pqr=True)
        fitter.go()
        res=fitter.get_result()

        return fitter

    def get_guess(self, pars, n=1):
        """
        Get a guess centered on the truth

        width is relative for T and counts
        """
        from ngmix.shape import g1g2_to_eta1eta2

        guess=numpy.zeros( (n, pars.size) )
        
        guess[:,0] = 0.001*srandu(n)
        guess[:,1] = 0.001*srandu(n)

        # convert pars to log
        eta1,eta2=g1g2_to_eta1eta2(pars[2],pars[3])
        logT=log10(pars[4])
        logF=log10(pars[5])

        guess[:,2] = eta1*(1.0 + 0.001*srandu(n))
        guess[:,3] = eta2*(1.0 + 0.001*srandu(n))
        guess[:,4] = logT*(1.0 + 0.01*srandu(n))
        guess[:,5] = logF*(1.0 + 0.01*srandu(n))

        return guess


    def get_pair_pars(self, **keys):
        """
        Get pair parameters
        """
        from numpy import pi
        import ngmix
        from ngmix.shape import eta1eta2_to_g1g2
        from ngmix.gexceptions import GMixRangeError

        # centers zero
        pars1=numpy.zeros(6)

        while True:
            try:
                logpars1=self.joint_prior.sample()

                eta1=logpars1[0]
                eta2=logpars1[1]
                logT=logpars1[2]
                logF=logpars1[3]
                
                # this might raise an exception
                g1,g2=eta1eta2_to_g1g2(eta1,eta2)
                T=10.0**(logT)
                F=10.0**(logF)

                pars1[2] = g1
                pars1[3] = g2
                pars1[4] = T
                pars1[5] = F
                break
            except GMixRangeError:
                pass

        print("logpars:",logpars1)

        pars2=pars1.copy()

        shape1 = ngmix.shape.Shape(pars1[2], pars1[3])
        shape2 = shape1.copy()
        shape2.rotate(pi/2.0)

        shear=self.shear
        shape1.shear(shear[0], shear[1])
        shape2.shear(shear[0], shear[1])

        pars1[2], pars1[3] = shape1.g1, shape1.g2
        pars2[2], pars2[3] = shape2.g1, shape2.g2

        cen_offset_arcsec=array( self.cen_prior.sample() )
        return pars1, pars2, cen_offset_arcsec

    def set_prior(self):
        """
        Set all the priors
        """
        import great3
        import ngmix


        self.pixel_scale = self.simc.get('pixel_scale',1.0)

        joint_prior_type = self.simc['joint_prior_type']
        self.joint_prior = \
            great3.joint_prior.make_joint_prior_simple(type=joint_prior_type)

        cen_sigma_arcsec=self.simc['cen_sigma']
        self.cen_prior=ngmix.priors.CenPrior(0.0,
                                             0.0,
                                             cen_sigma_arcsec,
                                             cen_sigma_arcsec)



def srandu(num=None):
    """
    Generate random numbers in the symmetric distribution [-1,1]
    """
    return 2*(numpy.random.random(num)-0.5)


def run_fmin(func, guess):
    import scipy.optimize
    (minvalx, fval, iterations, fcalls, warnflag) \
            = scipy.optimize.fmin(func,
                                  guess,
                                  full_output=True, 
                                  disp=False)
    if warnflag != 0:
        raise ValueError("failed to find min: warnflag %d" % warnflag)

    res={'flags':warnflag,
         'pars':minvalx,
         'nfev':fcalls}
    return res
def run_fmin_powell(func, guess):
    import scipy.optimize

    (minvalx, fval, direc, iterations, fcalls, warnflag) \
            = scipy.optimize.fmin_powell(func,
                                         guess,
                                         ftol=1.0e-6,
                                         xtol=1.0e-6,
                                         full_output=True, 
                                         disp=False)
    #if warnflag != 0:
    #    raise ValueError("failed to find min: warnflag %d" % warnflag)

    res={'flags':warnflag,
         'pars':minvalx,
         'fval':fval,
         'nfev':fcalls}
    return res

def test_simple_sampler():
    import ngmix
    import mcmc

    cen_dist    = ngmix.priors.Student2D(5.0, -3.0, 0.5, 0.5) 
    g_dist      = ngmix.priors.TruncatedStudentPolar(0.5, -0.8, 0.3, 0.3, 1.0)
    T_dist      =ngmix.priors.StudentPositive(16.0, 3.0)
    counts_dist =ngmix.priors.StudentPositive(100.0, 20.0)

    ss=SimpleSampler(cen_dist, g_dist, T_dist, counts_dist)
    samples = ss.sample(10000)

    mcmc.plot_results(samples)

class SimpleSampler(object):
    def __init__(self, cen_dist, g_dist, T_dist, counts_dist):
        """
        Each dist must have methods
            sample(n)
            get_lnprob_array(vals)
        """
        self.cen_dist=cen_dist
        self.g_dist=g_dist
        self.T_dist=T_dist
        self.counts_dist=counts_dist

        self.npars=6

    def sample(self, n):
        samples=numpy.zeros( (n, self.npars) )
        lnprob=numpy.zeros(n)

        cen1,cen2=self.cen_dist.sample(n)
        g1,g2 = self.g_dist.sample(n)
        T=self.T_dist.sample(n)
        counts=self.counts_dist.sample(n)

        samples[:, 0]=cen1
        samples[:, 1]=cen2
        samples[:, 2]=g1
        samples[:, 3]=g2
        samples[:, 4]=T
        samples[:, 5]=counts

        return samples

    def get_lnprob(self, pars):
        n=pars.shape[0]
        lnprob = numpy.zeros(n)

        lnprob += self.cen_dist.get_lnprob_array(pars[:,0],pars[:,1])
        lnprob += self.g_dist.get_lnprob_array(pars[:,2],pars[:,3])

        lnprob += self.T_dist.get_lnprob_array(pars[:,4])

        lnprob += self.counts_dist.get_lnprob_array(pars[:,5])

        return lnprob

    def get_prob(self, pars):
        lnprob = self.get_lnprob(pars)
        return numpy.exp( lnprob )


    def copy(self):
        return SimpleSampler(self.cen_dist,
                             self.g_dist,
                             self.T_dist,
                             self.counts_dist)
    def re_center(self, cen1, cen2, g_mean1, g_mean2, T_mean, counts_mean):
        d=self.cen_dist
        d.reset(cen1, cen2, d.sigma1, d.sigma2)

        d=self.g_dist
        d.reset(g_mean1, g_mean2, d.sigma1, d.sigma2, d.maxval)

        d=self.T_dist
        d.reset(T_mean, d.sigma)

        d=self.counts_dist
        d.reset(counts_mean, d.sigma)

def write_fits(filename, data):
    """
    Assume condor where cwd is scratch dir

    Write to cwd assuming scratch the move the file.

    The move may fail; it is retried a few times.
    """
    import fitsio
    import time

    output_file=os.path.abspath(filename)

    local_file=os.path.abspath( os.path.basename(output_file) )
    print("writing local file:",local_file)

    with fitsio.FITS(local_file,'rw',clobber=True) as fobj:
        fobj.write(data)

    if local_file==output_file:
        return True

    # remove if it exists
    try:
        os.remove(output_file)
    except:
        pass

    # try a few times
    print("moving to:",output_file)
    cmd='mv -v %s %s' % (local_file, output_file)
    for i in xrange(5):
        stat=os.system(cmd)
        if stat==0:
            print('success')
            success=True
            break
        else:
            print('error moving file, trying again')
            time.sleep(5)
            success=False

    return success

def seed_global_devrand():
    """
    Seed the "global" random number generator
    """
    seed = get_devrand_uint()
    numpy.random.seed(seed)

def get_random_state_devrand():
    """
    Seed the numpy random state from /dev/random
    """
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
