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
from ngmix.gexceptions import GMixRangeError
from ngmix.observation import Observation
from ngmix.gexceptions import GMixMaxIterEM

# region over which to render images and calculate likelihoods
NSIGMA_RENDER=5.0

# minutes
DEFAULT_CHECKPOINTS=[30,60,90,110]
#DEFAULT_CHECKPOINTS=[30,90]
#DEFAULT_CHECKPOINTS=[0,100]

from .fitters import TryAgainError
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
        self['verbose']=self.get('verbose',True)

        self._set_shear()
        self.nsub=self.simc['nsub']
        self.nsigma_render=self.simc.get('nsigma_render',NSIGMA_RENDER)

        self.check_pqr_shear()

        # this might be mean for correspond to a flux threshold, 
        # depending on the type of run
        self.s2n=s2n
        self.npairs=npairs

        self.create_model_set()

        self.fit_model=self['fit_model']


        if self.fit_model is not None:
            self.npars=ngmix.gmix.get_model_npars(self.fit_model)
            print("npars:",self.npars)

        self.make_plots=keys.get('make_plots',False)
        self.separate=False
        self.plot_base=keys.get('plot_base',None)

        self.setup_checkpoints(**keys)

        if self.data is None:
            self.make_struct()

        self.set_prior()
        self.set_dims_cen()
        self.make_psf()

        self.set_noise()

        if self['verbose']:
            pprint.pprint(self)
            pprint.pprint(self.simc)

    def _set_shear(self):
        if isinstance(self.simc['shear'], dict):
            shc=self.simc['shear']
            if shc['type']=='normal':
                self._shear_dist=ngmix.priors.RoundGauss2D(shc['mean'][0],
                                                           shc['mean'][1],
                                                           shc['sigma'][0],
                                                           shc['sigma'][1])
            else:
                raise ValueError("only lognormal for now")
            self['shear_is_constant'] = False
        else:
            self._shear_true=self.simc['shear']
            self['shear_is_constant'] = True

    def get_shear(self):
        """
        if shear is constant, return that, otherwise sample it
        """
        if self['shear_is_constant']:
            shear = self._shear_true
        else:
            shear = self._shear_dist.sample()

        return shear


    def create_model_set(self):
        """
        crae
        """
        from .sample import DiscreteSampler

        obj_model=self.simc['obj_model']
        if isinstance(obj_model,dict):
            sample_dict = obj_model
        else:
            sample_dict={obj_model: 100}

        self.model_sampler = DiscreteSampler(sample_dict)

        tmod = sample_dict.keys()[0]
        self.true_npars = ngmix.gmix.get_model_npars(tmod)

    def run_sim(self):
        """
        Run the simulation, fitting psf and all pairs
        """

        self.start_timer()

        self.index = 0
        i=self.index

        npairs=self.npairs
        for ipair in xrange(npairs):
            if self['verbose'] or ( (ipair % 100) == 0):
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
                        if self['verbose']:
                            print(str(err))

                self.copy_to_output(reslist[0], i)
                i += 1
                self.copy_to_output(reslist[1], i)
                i += 1

            self.set_elapsed_time()
            self.try_checkpoint()

        self.set_elapsed_time()

        print('time minutes:',self.tm_minutes)
        print('time per pair sec:',self.tm/npairs)
        print('time per image sec:',self.tm/(2*npairs))



    def set_config(self, sim_conf, run_conf):
        """
        Check and set the configurations
        """

        run_conf['use_logpars']=run_conf.get('use_logpars',False)

        if sim_conf['name'] != run_conf['sim']:
            err="sim name in run config '%s' doesn't match sim name '%s'"
            raise ValueError(err % (run_conf['sim'],sim_conf['name']))

        sim_conf['do_ring'] = sim_conf.get('do_ring',True)
        if not sim_conf['do_ring']:
            print("    not doing ring")

        self.simc=sim_conf
        self.conf=run_conf

    def check_pqr_shear(self):
        """
        If we are doing pqr need to check what shear we expand
        about
        """

        expand_shear_true=self.get('expand_shear_true',False)
        if expand_shear_true:
            self.shear_expand=array(self.simc['shear'])
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

        if self['verbose']:
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
            #res['pars_noshear'] = imd['pars_noshear']
            res['s2n_true'] = imd['s2n']
            if self['verbose']:
                self.print_res(res)

            if 'maxlike' in self['guess_type']:
                res['maxlike_res'] = self.maxlike_res

            self.fit_psf_flux(obs, res)
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

        pdict=fitter.make_plots(do_residual=True,
                                title=self.fit_model,
                                separate=self.separate,
                                weights=self._weights)

        if hasattr(self,'max_fitter'):
            resid_pname=self.plot_base+'-%06d-%s-gal-resid.png' % (self.ipair,key)
            print(resid_pname)
            rplt=self.max_fitter.plot_residuals()
            #pdict['resid'].write_img(width,height,resid_pname)
            rplt[0][0].write_img(width,height,resid_pname)
        
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
        if fitter_type in ['mcmc','pqrs-emcee']:
            fitter = self.fit_galaxy_mcmc(imdict)
        elif 'mh' in fitter_type:
            fitter = self.fit_galaxy_mh(imdict)
        elif fitter_type == 'lm':
            ntry=self['lm_ntry']
            fitter = self.fit_galaxy_lm(imdict,
                                        ntry=ntry)
            res=fitter.get_result()
            if res['flags']==0 and 'g' in res:
                check_g(res['g'])

        elif fitter_type=='lm-metacal':
            fitter=self.fit_galaxy_metacal(imdict, method='lm')
        elif fitter_type=='nm-metacal':
            fitter=self.fit_galaxy_metacal(imdict,method='nm')
        elif fitter_type=='em-metacal':
            fitter=self.fit_galaxy_metacal(imdict,method='em')
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
        g_prior=self.g_prior

        if self.g_prior_during:
            weights = None
            remove_prior=True
        else:
            weights = g_prior.get_prob_array2d(g[:,0], g[:,1])
            remove_prior=False

        # keep for later if we want to make plots
        self._weights=weights

        fitter.calc_result(weights=weights)

        # we are going to mutate the result dict owned by the fitter
        res=fitter.get_result()

        ls=ngmix.lensfit.LensfitSensitivity(g, g_prior,
                                            remove_prior=remove_prior)
        g_sens = ls.get_g_sens()
        g_mean = ls.get_g_mean()

        nuse = ls.get_nuse()

        pqrobj=ngmix.pqr.PQR(g, g_prior,
                             shear_expand=self.shear_expand,
                             remove_prior=remove_prior)


        P,Q,R = pqrobj.get_pqr()

        # this nuse should be the same for both lensfit and pqr
        res['nuse'] = nuse
        res['g_sens'] = g_sens
        res['P']=P
        res['Q']=Q
        res['R']=R


    def get_guess(self, imdict, n=1, guess_type=None):
        """
        Get a guess
        """

        if guess_type is None:
            guess_type=self['guess_type']

        if guess_type=='draw_truth':
            print('    * guessing randomized truth')
            guess=self.get_guess_from_pars(imdict['pars'], n=n)
        elif guess_type=='draw_pdf':
            guess=self.get_guess_draw_pdf(n=n)
        elif guess_type=='draw_maxlike':
            guess,perr=self.get_guess_draw_maxlike(imdict, n=n)
        else:
            raise ValueError("bad guess type: '%s'" % guess_type)
        
        if self['use_logpars'] and guess_type in ['draw_truth','draw_pdf']:
            if len(guess.shape)==1:
                guess[4:] = log(guess[4:])
            else:
                guess[:, 4:] = log(guess[:, 4:])
        return guess


    def run_simple_mcmc_fitter(self, imdict):
        """
        simple gauss,exp,dev
        """
        from ngmix.fitting import MCMCSimple

        epars = self['emcee_pars']

        guess=self.get_guess(imdict, n=epars['nwalkers'])

        obs=imdict['obs']

        fitter=MCMCSimple(obs,
                          self.fit_model,

                          prior=self.search_prior, # no prior during

                          nwalkers=epars['nwalkers'],
                          mca_a=epars['a'],
                          random_state=self.random_state,
                          use_logpars=self['use_logpars'])

        pos=fitter.run_mcmc(guess,epars['burnin'])
        pos=fitter.run_mcmc(pos,epars['nstep'], thin=epars['thin'])

        if 'min_arate' in epars:
            arate=fitter.get_arate()
            if arate < epars['min_arate']:
                print("        low arate",arate,"running some more")
                pos=fitter.run_mcmc(pos,epars['burnin']*2)
                pos=fitter.run_mcmc(pos,epars['nstep'], thin=epars['thin'])


        return fitter

    def run_simple_mh_fitter(self, imdict):
        """
        simple gauss,exp,dev
        """
        from ngmix.fitting import MHSimple, MHTempSimple

        mess="for mh guess should be maxlike"
        assert (self['guess_type']=="draw_maxlike"),mess

        guess,perr=self.get_guess(imdict, n=1, guess_type='draw_maxlike')


        obs=imdict['obs']

        temp=self.get('temp',None)
        if temp is not None:
            print("run_simple_mh_fitter: doing temperature:",temp)
            # note modified step sizes for temperature
            step_sizes = 0.5*perr*sqrt(temp)
            fitter=MHTempSimple(obs,
                                self.fit_model,
                                step_sizes,

                                temp=temp,
                                use_logpars=self['use_logpars'],

                                prior=self.search_prior, 
                                random_state=self.random_state)

        else:
            step_sizes = 0.5*perr
            fitter=MHSimple(obs,
                            self.fit_model,
                            step_sizes,

                            prior=self.search_prior, 
                            random_state=self.random_state)

        pos=fitter.run_mcmc(guess,self['burnin'])
        pos=fitter.run_mcmc(pos,self['nstep'])

        return fitter


    def get_guess_draw_pdf(self, n=1):
        """
        Get a guess drawn from the priors

        assume simple for now
        """
        if self['verbose']:
            print("drawing guess from priors")
        guess=self.pdf.sample(n)


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

        max_guess=self.get_guess(imdict, n=1, guess_type='draw_pdf')

        fitter=self.fit_galaxy_max(imdict, max_guess)

        res=fitter.get_result()
        self.maxlike_res=res

        pars=res['pars']
        perr=res.get('pars_err',None)

        print("    nfev:",res['nfev'])
        print_pars(pars, front='    max pars: ')

        if perr is not None:
            print_pars(perr, front='    max perr: ')

            # try diagonal for now

            width = perr.copy()
            width[0] = width[0].clip(min=0.001,max=1.0)
            width[1] = width[1].clip(min=0.001,max=1.0)
            width[2] = width[2].clip(min=0.001,max=0.5)
            width[3] = width[3].clip(min=0.001,max=0.5)
            width[4] = width[4].clip(min=0.01,max=100.0)
            width[5] = width[5].clip(min=0.1,max=100.0)

            tcov = diag(width**2)
            guess=self.get_guess_from_cov(pars,
                                          tcov,
                                          self.search_prior,
                                          nrand=n)
        else:
            print("    no cov")
            # now scatter it around a bit
            # if logpars these are also logpars
            guess=self.get_guess_from_pars(pars, n=n)


        return guess, perr

    def get_guess_from_pars(self, pars, n=None, width=None):
        """
        Get a guess centered on the input linear pars

        uniform sampling
        width is fractional for T and flux
        """

        if n is None:
            is_scalar=True
            n=1
        else:
            is_scalar=False
        npars=pars.size

        if width is None:
            width = pars*0 + 0.01
        else:
            if len(width) != len(pars):
                raise ValueError("width not same size as pars")

        guess=numpy.zeros( (n, pars.size) )

        guess[:,0] = pars[0] + width[0]*srandu(n)
        guess[:,1] = pars[1] + width[1]*srandu(n)

        guess_shape=self.get_shape_guess(pars[2],pars[3],n,width=width[2:2+2])
        guess[:,2]=guess_shape[:,0]
        guess[:,3]=guess_shape[:,1]

        for i in xrange(4,npars):
            guess[:,i] = pars[i]*( 1.0 + width[i]*srandu(n) )

        if is_scalar:
            guess=guess[0,:]

        return guess

    def get_guess_from_cov(self, pars, cov, prior, nrand=None):
        """
        Get a guess centered on the input pars

        you should clip yourself
        """
        from numpy.random import multivariate_normal
        
        if nrand is None:
            nrand=1
            is_scalar=True
        else:
            is_scalar=False

        guess=numpy.zeros( (nrand, pars.size) )

        for i in xrange(nrand):

            while True:
                try:
                    rpars=multivariate_normal(pars, cov)
                    lnp = prior.get_lnprob_scalar(rpars)
                    break
                except GMixRangeError:
                    pass

            guess[i,:] = rpars

        if is_scalar:
            guess=guess[0,:]

        return guess


    def get_shape_guess(self, g1, g2, n, width=[0.01,0.01], rsampler=None):
        """
        Get guess, making sure in range
        """

        if rsampler is None:
            rsampler=srandu

        guess=numpy.zeros( (n, 2) )
        shape=ngmix.Shape(g1, g2)

        for i in xrange(n):

            while True:
                try:
                    g1_offset = width[0]*rsampler()
                    g2_offset = width[1]*rsampler()
                    shape_new=shape.copy()
                    shape_new.shear(g1_offset, g2_offset)
                    break
                except GMixRangeError:
                    pass

            guess[i,0] = shape_new.g1
            guess[i,1] = shape_new.g2

        return guess


    def fit_galaxy_max(self, imdict, guess):
        """
        Fit the model to the galaxy

        we send some keywords so behavior can be different if this is a guess
        getter

        """
        from ngmix.fitting import MaxSimple

        obs=imdict['obs']
        fitter=MaxSimple(obs, self.fit_model, prior=self.search_prior,
                         use_logpars=self['use_logpars'])
        fitter.run_max(guess, **self['nm_pars'])
        res=fitter.get_result()
        res['ntry']=1
        return fitter



    def fit_galaxy_lm(self, imdict, ntry=1, guess=None):
        """
        Fit the model to the galaxy

        we send some keywords so behavior can be different if this is a guess
        getter

        """
        from ngmix.fitting import LMSimple

        obs=imdict['obs']
        fitter=LMSimple(obs,
                        self.fit_model,
                        prior=self.search_prior,
                        lm_pars=self['lm_pars'],
                        use_logpars=self['use_logpars'])

        for i in xrange(ntry):
            if i > 0 or guess is None:
                guess=self.get_lm_guess(imdict)

            fitter.run_lm(guess)
            res=fitter.get_result()
            if res['flags']==0:
                break

        res['ntry']=i+1

        #if 'pars' in res:
        #    print_pars(res['pars'], front='      metacal: ')

        if res['flags']==0:
            print("        replacing cov")
            # try to replace the covariance
            h=1.0e-3
            m=5.0
            fitter.calc_cov(h, m)
            if res['flags'] != 0:
                print("        replacement failed")
                res['flags']=0

        return fitter


    def fit_galaxy_nm(self, imdict, ntry=1, guess=None):
        """
        Fit the model to the galaxy

        we send some keywords so behavior can be different if this is a guess
        getter

        """
        from ngmix.fitting import MaxSimple

        nm_pars=self['nm_pars']

        obs=imdict['obs']
        fitter=MaxSimple(obs, self.fit_model,
                         method='nm',
                         prior=self.search_prior,
                         use_logpars=self['use_logpars'])

        for i in xrange(ntry):
            if i > 0 or guess is None:
                guess=self.get_lm_guess(imdict)

            fitter.go(guess, **nm_pars)
            res=fitter.get_result()
            if res['flags']==0:
                break

        if 'pars' in res:
            print_pars(res['pars'], front='      metacal: ')
        res['ntry']=i+1
        return fitter

    def fit_galaxy_em(self, imdict, ntry=1):
        """
        Fit the model to the galaxy

        we send some keywords so behavior can be different if this is a guess
        getter

        """
        from ngmix.em import GMixEM, prep_image

        obs=imdict['obs']

        image=obs.image
        imsky,sky=prep_image(image)

        sky_obs = Observation(imsky, jacobian=obs.jacobian)

        fitter=GMixEM(sky_obs)

        tol=1.0e-6
        maxiter=30000

        for i in xrange(ntry):
            guess=self.get_psf_guess()
            try:
                fitter.go(guess, sky, tol=tol,maxiter=maxiter)
                break
            except GMixMaxIterEM as e:
                print(str(e))
                print('re-trying em gal')

        res=fitter.get_result()
        if res['flags'] != 0:
            raise TryAgainError("em gal failed")


        gm=fitter.get_gmix()
        pars=gm.get_full_pars()

        T=pars[3]+pars[5]
        e1 = (pars[5]-pars[3])/T
        e2 = 2*pars[4]/T
        g1,g2=e1,e2
        #try:
        #    g1,g2=ngmix.shape.e1e2_to_g1g2(e1,e2)
        #except GMixRangeError:
        #    raise TryAgainError("bad e1,e2")

        tpars=array([pars[1],pars[2],
                     e1,e2,
                     T,
                     1.0])

        res['g']=array([g1,g2])
        res['g_cov']=diag([9999.0,9999.0])
        res['pars']=tpars
        res['pars_err']=pars*0 + 9999.0
        res['pars_cov']=diag(pars*0 + 9999)
        res['s2n_w']=-9999.0
        res['nfev']=res['numiter']
        res['ntry']=i+1

        return fitter

    def fit_galaxy_metacal(self, imdict, method='lm'):
        """
        Fit the model to the galaxy, and also calculate a sensitivity
        """
        from ngmix.fitting import LMSimple

        mpars=self['metacal_pars']

        ntry=self['lm_ntry']
        if method=='lm':
            fmethod=self.fit_galaxy_lm
        elif method=='em':
            fmethod=self.fit_galaxy_em
        else:
            fmethod=self.fit_galaxy_nm

        fitter=fmethod(imdict, ntry=ntry)

        res=fitter.get_result()
        if res['flags'] != 0:
            raise TryAgainError("bad lm measurements")

        # might raise exception
        check_g(res['g'])

        imdict_lo, imdict_hi = self.get_metacal_imdicts(imdict, fitter)

        fitter_lo=fmethod(imdict_lo, ntry=ntry)
        fitter_hi=fmethod(imdict_hi, ntry=ntry)

        res_lo=fitter_lo.get_result()
        res_hi=fitter_hi.get_result()

        if res_lo['flags'] != 0 or res_hi['flags'] != 0:
            raise TryAgainError("bad metacal measurements")

        check_g(res_lo['g'])
        check_g(res_hi['g'])

        pars_lo=res_lo['pars']
        pars_hi=res_hi['pars']

        # assuming same for now
        h=self['metacal_pars']['step']
        g_sens = array( [(pars_hi[2]-pars_lo[2])/(2.*h)]*2 )

        print("        sensitivity: %g %g" % (g_sens[0],g_sens[1]))
        res['g_sens'] = g_sens

        return fitter

    def get_metacal_imdicts(self, imdict, fitter):

        mpars=self['metacal_pars']
        h=mpars['step']
        pars2use=mpars['pars2use']

        obs_orig=imdict['obs']

        psf_gmix, pars, model, nsub=self.get_metacal_input(imdict,fitter,pars2use)

        # only doing shear in g1 direction since that is what
        # we usually have for sims. Can adapt later to real data
        pars_lo=make_sheared_pars(pars, -h, 0.0)
        pars_hi=make_sheared_pars(pars, +h, 0.0)

        im=obs_orig.image

        noise2use=mpars['noise2use']
        if noise2use=='same':
            noise_image = imdict['nim']
        elif noise2use in ['same-sim','independent']:
            noise_image = self.skysig*randn(im.size).reshape(im.shape)

        obs_lo = make_metacal_obs(psf_gmix,
                                  pars_lo, 
                                  model,
                                  noise_image,
                                  obs_orig,
                                  nsub)

        if noise2use=='independent':
            noise_image = self.skysig*randn(im.size).reshape(im.shape)

        obs_hi = make_metacal_obs(psf_gmix,
                                  pars_hi, 
                                  model,
                                  noise_image,
                                  obs_orig,
                                  nsub)

        imdict_lo = {'obs':obs_lo,'pars':pars_lo}
        imdict_hi = {'obs':obs_hi,'pars':pars_hi}

        return imdict_lo, imdict_hi

    def get_metacal_input(self, imdict, fitter, pars2use):

        if pars2use=='fit':
            # render a sheared version of our fit, note nsub=1
            psf=imdict['obs'].psf.gmix.copy()
            pars=fitter.get_result()['pars'].copy()
            model=self.fit_model
            nsub=1
        else:
            # Use the original model but sheared differently.
            # closer to what Eric suggests with data manipulations

            psf=self.psf_gmix_true.copy()
            pars=imdict['pars'].copy()
            model=self.simc['obj_model']
            nsub=self.nsub

        return psf, pars, model, nsub


    def get_lm_guess(self, imdict):
        """
        Get a guess for the LM fitter
        """
        guess_type=self['guess_type']
        if guess_type=='truth':
            guess=imdict['pars']
        elif guess_type=='truth_random':
            guess=self.get_guess_from_pars(imdict['pars'])
        elif guess_type=='draw_pdf':
            guess=self.get_guess_draw_pdf(n=1)
        else:
            raise ValueError("bad guess type: '%s'" % guess_type)
        return guess

    def print_res(self,res):
        """
        print some stats
        """

        if 'arate' in res:
            mess="    arate: %.2f tau: %.1f s2n_w: %.1f nuse: %d"
            mess=mess % (res['arate'],res['tau'],res['s2n_w'],res['nuse'])
            print(mess)
        elif 'nfev' in res:
            print('    try:',res['ntry'],'nfev:',res['nfev'],'s2n_w:',res['s2n_w'])

        print_pars(res['pars'],front='    pars: ')
        print_pars(res['pars_err'],front='    perr: ')
        print_pars(res['pars_true'], front='    true: ')

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

        image=obs.image
        imsky,sky=ngmix.em.prep_image(image)

        sky_obs = Observation(imsky, jacobian=obs.jacobian)

        em=ngmix.em.GMixEM(sky_obs)

        tol=self.get('psf_tol',1.0e-6)
        maxiter=self.get('psf_maxiter',30000)

        while True:
            guess=self.get_psf_guess()
            try:
                em.go(guess, sky, tol=tol,maxiter=maxiter)
                break
            except GMixMaxIterEM as e:
                print(str(e))
                print('re-trying psf')

        psf_gmix_fit=em.get_gmix()

        if self['verbose']:
            print('psf fit:')
            print(psf_gmix_fit)

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
            
            guess=ngmix.gmix.GMix(pars=pars)

        return guess

    def fit_psf_flux(self, obs, res):
        """
        use psf as a template, measure flux (linear)
        """

        if True:
            cen=res['pars'][0:2]
            fitter=ngmix.fitting.TemplateFluxFitter(obs, do_psf=True, cen=cen)
        elif False:
            fitter=ngmix.fitting.TemplateFluxFitter(obs, do_psf=True)
        else:
            pars=res['pars'].copy()
            pars[4] = exp(pars[4])
            pars[5] = exp(pars[5])
            gm0=ngmix.GMixModel(res['pars'], self['fit_model'])
            gm=gm0.convolve(obs.psf.gmix)
            obs.set_gmix(gm)

            fitter=ngmix.fitting.TemplateFluxFitter(obs)

        fitter.go()

        pres=fitter.get_result()


        res['psf_flux'] = pres['flux']
        res['psf_flux_err'] = pres['flux_err']
        res['psf_flux_s2n'] = pres['flux']/pres['flux_err']

        print("    psf flux: %.3f +/- %.3f" % (pres['flux'],pres['flux_err']))




    def set_dims_cen(self):
        """
        When using constant dims
        """
        from math import ceil

        if 'box_size' in self.simc:
            self.dims_pix=array( [self.simc['box_size']]*2, dtype='i4')
        else:
            samples = self.pdf.sample(10000)
            T = samples[:,4]

            Tstd = T.std()

            Tmax = self.simc['obj_T_mean'] + 3.0*Tstd

            sigma_pix=numpy.sqrt(Tmax/2.)/self.pixel_scale
            dim=int(ceil(2.*sigma_pix*self.nsigma_render))
            self.dims_pix = array( [dim,dim], dtype='i8')

        self.cen_pix = array( [(self.dims_pix[0]-1.)/2.]*2 )

        print("dims:",self.dims_pix)
        print("cen: ",self.cen_pix)

    def set_prior(self):
        """
        Set all the priors
        """
        import ngmix
        from ngmix.joint_prior import PriorSimpleSep

        self.pixel_scale = self.simc.get('pixel_scale',1.0)

        # we may over-ride below
        self.s2n_for_noise = self.s2n

        simc=self.simc

        # prior separate for all pars
        cen_sigma_arcsec=simc['cen_sigma']
        cen_pdf=ngmix.priors.CenPrior(0.0,
                                      0.0,
                                      cen_sigma_arcsec,
                                      cen_sigma_arcsec)

        gtype=simc['g_prior_type']
        if gtype=="ba":
            g_pdf_sigma=simc['g_prior_sigma']
            g_pdf=ngmix.priors.GPriorBA(g_pdf_sigma)
        elif gtype=="cosmos":
            g_pdf=ngmix.priors.make_gprior_cosmos_sersic(type='erf')
        elif gtype=='great-des':
            g_pdf= ngmix.priors.GPriorGreatDES(pars=simc['g_prior_pars'],
                                               gmax=1.0)

        else:
            raise ValueError("only g prior 'ba' for now")

        g_pdf_flat=ngmix.priors.ZDisk2D(1.0)

        # T and scatter in linear space
        T            = simc['obj_T_mean']
        T_sigma      = simc['obj_T_sigma_frac']*T
        counts       = simc['obj_counts_mean']
        counts_sigma = simc['obj_counts_sigma_frac']*counts


        T_pdf=ngmix.priors.LogNormal(T, T_sigma)
        counts_pdf=ngmix.priors.LogNormal(counts, counts_sigma)

        # for drawing parameters, and after exploration to grab g_pdf and
        # calculate pqr etc.
        self.pdf = PriorSimpleSep(cen_pdf,
                                  g_pdf,
                                  T_pdf,
                                  counts_pdf)


        gppars = self['g_prior_pars']
        if gppars['type']=='ba':
            self.g_prior = ngmix.priors.GPriorBA(gppars['sigma'])
        elif gppars['type']=='cosmos':
            self.g_prior=ngmix.priors.make_gprior_cosmos_sersic(type='erf')
        elif gppars['type']=='great-des':
            self.g_prior = ngmix.priors.GPriorGreatDES(pars=gppars['pars'],
                                                       gmax=1.0)

        else:
            raise ValueError("implement other")

        self.g_prior_during = self.get('g_prior_during',False)
        if self.g_prior_during:
            print("    Using full g prior for fits")
            fit_g_prior = self.g_prior
        else:
            print("    Using flat g prior for fits")
            fit_g_prior = g_pdf_flat

        print("using input search prior")
        spars=self['search_prior']

        if spars['T_prior_type']=="truth":
            print("using true T pdf for prior")
            if self['use_logpars']:
                print("    converting to log")
                logT_mean, logT_sigma=ngmix.priors.lognorm_convert(T,T_sigma)
                fit_T_prior = ngmix.priors.Normal(logT_mean, logT_sigma)
            else:
                fit_T_prior = T_pdf
        else:
            T_prior_pars = spars['T_prior_pars']
            fit_T_prior=ngmix.priors.TwoSidedErf(*T_prior_pars)

        if spars['counts_prior_type']=="truth":
            print("using true counts pdf for prior")
            if self['use_logpars']:
                print("    converting to log")
                logc_mean, logc_sigma=ngmix.priors.lognorm_convert(counts,
                                                                   counts_sigma)
                fit_counts_prior = ngmix.priors.Normal(logc_mean, logc_sigma)
            else:
                fit_counts_prior = counts_pdf
            #fit_counts_prior = counts_pdf
        else:
            counts_prior_pars = spars['counts_prior_pars']
            fit_counts_prior=ngmix.priors.TwoSidedErf(*counts_prior_pars)

        if spars['cen_prior_type']=="truth":
            fit_cen_prior=cen_pdf
        else:
            fit_cen_sigma_arcsec=spars['cen_sigma']
            fit_cen_prior=ngmix.priors.CenPrior(0.0,
                                                0.0,
                                                fit_cen_sigma_arcsec,
                                                fit_cen_sigma_arcsec)


        self.search_prior = PriorSimpleSep(fit_cen_prior,
                                           fit_g_prior,
                                           fit_T_prior,
                                           fit_counts_prior)

    
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

            imn, noise_im =self.add_noise(im)
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
        im1,nim1=self.add_noise(im1_0)
        im2,nim2=self.add_noise(im2_0)

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
        imdict['im1']['nim'] = nim1
        imdict['im2']['nim'] = nim2

        im1_s2n = self.get_model_s2n(im1_0)
        im2_s2n = self.get_model_s2n(im2_0)

        imdict['im1']['s2n'] = im1_s2n
        imdict['im2']['s2n'] = im2_s2n

        return imdict

    def add_noise(self, im):
        """
        Add gaussian random noise
        """
        nim = self.skysig*randn(im.size).reshape(im.shape)
        return im + nim, nim

    def get_image_pair(self, random=True):
        """
        get a model image

        If random is True, use draw random values from the priors.
        Otherwise use the mean of the priors

        """
        import ngmix
        from ngmix.observation import Observation

        pars1, pars2, p1_noshear, p2_noshear = self.get_pair_pars(random=random)

        obj_model = self.model_sampler()
        self.current_model_true = obj_model
        print("sim model:",obj_model)
        gm1_pre=ngmix.gmix.GMixModel(pars1, obj_model)
        gm2_pre=ngmix.gmix.GMixModel(pars2, obj_model)

        gm1  = gm1_pre.convolve(self.psf_gmix_true)
        gm2  = gm2_pre.convolve(self.psf_gmix_true)

        # in case not doing a ring
        T = max(gm1.get_T(), gm2.get_T())

        dims_pix, cen0_pix = self.get_dims_cen_pergal(T + self.simc['psf_T'])

        # conversion between pixels and sky in arcsec
        self.jacobian=ngmix.jacobian.Jacobian(cen0_pix[0],
                                              cen0_pix[1],
                                              self.pixel_scale,
                                              0.0,
                                              0.0,
                                              self.pixel_scale)

        psf_cen=pars1[0:0+2].copy()
        psf_image=self.get_psf_image(dims_pix, psf_cen)

        psf_obs = Observation(psf_image, jacobian=self.jacobian)

        nsub = self.nsub
        im1=gm1.make_image(dims_pix, nsub=nsub, jacobian=self.jacobian)
        im2=gm2.make_image(dims_pix, nsub=nsub, jacobian=self.jacobian)

        obs1 = Observation(im1, psf=psf_obs, jacobian=self.jacobian)
        obs2 = Observation(im2, psf=psf_obs, jacobian=self.jacobian)
        
        out={'im1':{'pars':pars1,
                    'pars_noshear':p1_noshear,
                    'gm_pre':gm1_pre,
                    'gm':gm1,
                    'obs0':obs1},
             'im2':{'pars':pars2,
                    'pars_noshear':p2_noshear,
                    'gm_pre':gm2_pre,
                    'gm':gm2,
                    'obs0':obs2}}
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

            if self.simc['do_ring']:
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

        shear=self.get_shear()
        #print("    shear:",shear)
        shape1.shear(shear[0], shear[1])
        shape2.shear(shear[0], shear[1])

        pars1[2]=shape1.g1
        pars1[3]=shape1.g2
        pars2[2]=shape2.g1
        pars2[3]=shape2.g2

        return pars1, pars2, pars1_noshear, pars2_noshear

    def get_dims_cen(self):
        """
        For fixed size boxes
        """

        return self.dims_pix.copy(), self.cen_pix.copy()

    def get_dims_cen_pergal(self, T):
        """
        Based on T, get the required dimensions and a center

        Should send T+Tpsf
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

        d['model_true'][i] = self.current_model_true


        d['s2n_true'][i] = res['s2n_true']
        d['pars_true'][i,:] = res['pars_true']
        #d['pars_noshear'][i,:] = res['pars_noshear']

        d['pars'][i,:] = res['pars']
        d['pcov'][i,:,:] = res['pars_cov']

        d['psf_flux'][i] = res['psf_flux']
        d['psf_flux_err'][i] = res['psf_flux_err']
        d['psf_flux_s2n'][i] = res['psf_flux_s2n']

        d['g'][i,:] = res['g']
        d['g_cov'][i,:,:] = res['g_cov']

        d['s2n_w'][i] = res['s2n_w']

        if 'g_sens' in res:
            d['g_sens'][i,:] = res['g_sens']

        if 'arate' in res:
            d['arate'][i] = res['arate']
            d['tau'][i] = res['tau']

        if 'P' in res:
            d['P'][i] = res['P']
            d['Q'][i,:] = res['Q'][:]
            d['R'][i,:,:] = res['R'][:,:]

            d['nuse'][i] = res['nuse']

        if 'nfev' in res:
            d['nfev'][i] = res['nfev']
            # set outside of fitter
            d['ntry'][i] = res['ntry']

        if 'maxlike' in self['guess_type']:
            mres=res['maxlike_res']
            d['max_pars'][i] = mres['pars']
            if 'chi2per' in mres:
                d['max_chi2per'][i] = mres['chi2per']

            if 'pars_cov' in mres:
                d['max_pcov'][i] = mres['pars_cov']

    def get_dtype(self):
        """
        get the dtype for the output struct
        """
        npars=self.npars

        true_npars = ngmix.gmix.get_model_npars(tmod)

        dt=[('processed','i2'),
            ('model_true','S3'),
            ('s2n_true','f8'),
            ('pars_true','f8',true_npars),
            ('pars','f8',npars),
            ('pcov','f8',(npars,npars)),
            ('g','f8',2),
            ('g_cov','f8',(2,2)),
            ('psf_flux','f8'),
            ('psf_flux_err','f8'),
            ('psf_flux_s2n','f8'),
            ('s2n_w','f8'),
            ('arate','f8'),
            ('tau','f8'),
           ]

        fitter=self['fitter']
        if 'lm' in fitter or 'nm' in fitter or 'em' in fitter:
            dt += [('nfev','i4'),
                   ('ntry','i4')]
            if 'metacal' in fitter:
                dt.append( ('g_sens','f8',2) )
        else:
            dt += [('g_sens','f8',2),
                   ('P','f8'),
                   ('Q','f8',2),
                   ('R','f8',(2,2)),
                   ('nuse','i4')]

        if 'maxlike' in self['guess_type']:
            dt += [('max_pars','f8',npars),
                   ('max_pcov','f8',(npars,npars)),
                   ('max_chi2per','f4')]

        return dt

    def make_struct(self):
        """
        Make the output array
        """

        dt=self.get_dtype()
        self.data=numpy.zeros(self.npairs*2, dtype=dt)

        if 'maxlike' in self['guess_type']:
            self.data['max_pcov'] = 9999.e9


class NGMixSimPQRSMCMC(NGMixSim):
    def _add_mcmc_stats(self, fitter):
        """
        Calculate some stats

        The result dict internal to the fitter is modified to include
        g_sens and P,Q,R

        could we add this as a stand-alone function to ngmix.fitting?
        """

        temperature_weights = fitter.get_weights()
        mess="pqrs doesn't support temperature weights right now"
        assert temperature_weights==None,mess

        trials = fitter.get_trials()
        g=trials[:,2:2+2]

        # this is the full prior
        g_prior=self.g_prior
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

        pqrs_obj=ngmix.pqr.PQRS(g, g_prior)

        P,Q,R,S = pqrs_obj.get_pqrs()

        # this nuse should be the same for both lensfit and pqr
        res['nuse'] = nuse
        res['g_sens'] = g_sens
        res['P']=P
        res['Q']=Q
        res['R']=R
        res['S']=S

    def copy_to_output(self, res, i):
        """
        Copy results into the output
        """
        super(NGMixSimPQRS,self).copy_to_output(res, i)
        self.data['S'][i,:,:,:] = res['S'][:,:,:]

    def get_dtype(self):
        """
        Augment the dtype of paranet with S
        """
        dt=super(NGMixSimPQRS,self).get_dtype()
        dt += [('S','f8',(2,2,2))]

        return dt

class NGMixSimPQRSLM(NGMixSim):
    def fit_galaxy(self, imdict):
        fitter = super(NGMixSimPQRSLM,self).fit_galaxy_lm(imdict,
                                                          ntry=self['lm_ntry'])

        res=fitter.get_result()
        if res['flags'] == 0:
            g=res['pars'][2:2+2]

            # this is the full prior
            g_prior=self.g_prior

            pqrs_obj=ngmix.pqr.PQRS(g, g_prior)

            P,Q,R,S = pqrs_obj.get_pqrs()

            # this nuse should be the same for both lensfit and pqr
            res['nuse'] = nuse
            res['g_sens'] = g_sens
            res['P']=P
            res['Q']=Q
            res['R']=R
            res['S']=S

        return fitter

    def copy_to_output(self, res, i):
        """
        Copy results into the output
        """
        super(NGMixSimPQRS,self).copy_to_output(res, i)
        self.data['S'][i,:,:,:] = res['S'][:,:,:]

    def get_dtype(self):
        """
        Augment the dtype of paranet with S
        """
        dt=super(NGMixSimPQRS,self).get_dtype()
        dt += [('S','f8',(2,2,2))]

        return dt



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
                                      use_logpars=self['use_logpars'],

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

        shear=self.get_shear()
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

                                     use_logpars=self['use_logpars'],

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

        pars1[4] = sample[0]
        pars1[5] = sample[1]

        shape1=ngmix.shape.Shape(g1a[0],g2a[0])

        shape2=shape1.copy()
        shape2.rotate(pi/2.0)

        shear=self.get_shear()

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

        guess[:,4] = pars[4]*(1.0 + width[4]*srandu(n))
        guess[:,5] = pars[5]*(1.0 + width[5]*srandu(n))


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

                               use_logpars=self['use_logpars'],

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

        shear=self.get_shear()
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



class NGMixSimISample(NGMixSim):
    def __init__(self, *args, **kw):
        super(NGMixSimISample,self).__init__(*args, **kw)
        
        #assert self.g_prior_during==True,"g prior during is required"

        ipars=self['isample_pars']
        ipars['min_err'] = array(ipars['min_err'])
        ipars['max_err'] = array(ipars['max_err'])

    def fit_galaxy(self, imdict):
        """
        do max and sample truncated gausian in g1,g2 for
        lensfit
        """


        ipars=self['isample_pars']

        max_fitter = self.run_max_fitter(imdict, self.fit_model)
        self.max_fitter=max_fitter

        use_fitter = max_fitter
        niter=len(ipars['nsample'])
        for i,nsample in enumerate(ipars['nsample']):

            sampler=self._make_sampler(use_fitter)
            sampler.make_samples(nsample)

            sampler.set_iweights(max_fitter.calc_lnprob)
            if self.g_prior_during:
                weights=None
            else:
                samples=sampler.get_samples()
                weights=self.g_prior.get_prob_array2d(samples[:,2],samples[:,3])
            sampler.calc_result(weights=weights)

            tres=sampler.get_result()

            print("    eff iter %d: %.2f" % (i,tres['efficiency']))
            use_fitter = sampler

        self._add_mcmc_stats(sampler)
        self._sampler=sampler

        res=sampler.get_result()
        res['maxlike_res']=max_fitter.get_result()
        return sampler

    def _add_mcmc_stats(self, sampler):
        """
        Calculate some stats

        The result dict internal to the sampler is modified to include
        g_sens and P,Q,R

        call calc_result before calling this method
        """

        # this is the full prior
        g_prior=self.g_prior

        iweights = sampler.get_iweights()
        samples = sampler.get_samples()
        g_vals=samples[:,2:2+2]

        if self.g_prior_during:
            remove_prior=True
        else:
            remove_prior=False

        # keep for later if we want to make plots
        self._weights=iweights

        # we are going to mutate the result dict owned by the sampler
        res=sampler.get_result()

        maxres = self.maxlike_res
        res['s2n_w'] = maxres['s2n_w']

        ls=ngmix.lensfit.LensfitSensitivity(g_vals,
                                            g_prior,
                                            weights=iweights,
                                            remove_prior=remove_prior)
        g_sens = ls.get_g_sens()
        g_mean = ls.get_g_mean()

        #print("        res g:",res['g'])
        #print("        lf g: ",g_mean)
        res['g_sens'] = g_sens
        res['nuse'] = ls.get_nuse()

        # not able to use extra weights yet
        pqrobj=ngmix.pqr.PQR(g_vals, g_prior,
                             shear_expand=self.shear_expand,
                             weights=iweights,
                             remove_prior=remove_prior)


        P,Q,R = pqrobj.get_pqr()
        res['P']=P
        res['Q']=Q
        res['R']=R
        

    def _make_sampler(self, fitter):
        from ngmix.fitting import ISampler
        from numpy.linalg import LinAlgError

        ipars=self['isample_pars']

        res=fitter.get_result()
        icov = res['pars_cov']*ipars['ifactor']**2

        try:
            sampler=ISampler(res['pars'],
                             icov,
                             ipars['df'],
                             min_err=ipars['min_err'],
                             max_err=ipars['max_err'])
        except LinAlgError:
            raise TryAgainError("bad cov")

        return sampler

    def run_max_fitter(self, imdict, model):
        """
        get the covariance matrix based on a max like fit

        raises if not found
        """

        ipars=self['isample_pars']

        if ipars['max_fitter']=="nm":
            fitmethod=self.fit_galaxy_nm
        else:
            fitmethod=self.fit_galaxy_lm

        for itry in xrange(1,ipars['max_tries']+1):
            max_guess=self.get_guess(imdict, guess_type='draw_truth')

            max_guess[4:] = log(max_guess[4:])

            fitter=fitmethod(imdict, max_guess, model)

            res=fitter.get_result()
            if res['flags']==0:
                break

        pars=res['pars']
        perr=res.get('pars_err',None)
        pars_cov=res.get('pars_cov',None)

        mess="    tries: %d  nfev: %d"
        mess=mess % (itry,res['nfev'])
        print(mess)
        print_pars(pars, front='    max pars: ')

        if res['flags'] != 0:
            if perr is None:
                print("    no cov")
            else:
                print("    did not converge")
            raise TryAgainError("bad max fit")

        print_pars(perr, front='    max perr: ')

        self.maxlike_res=res

        return fitter

    def fit_galaxy_nm(self, imdict, guess, model):
        """
        Fit the model to the galaxy

        we send some keywords so behavior can be different if this is a guess
        getter

        """
        from ngmix.fitting import MaxSimple

        obs=imdict['obs']
        fitter=MaxSimple(obs, model,
                         prior=self.search_prior,
                         use_logpars=self['use_logpars'])

        fitter.run_max(guess, **self['nm_pars'])
        return fitter

    def fit_galaxy_lm(self, imdict, guess, model):
        """
        Fit the model to the galaxy

        we send some keywords so behavior can be different if this is a guess
        getter

        """
        from ngmix.fitting import LMSimple

        obs=imdict['obs']
        fitter=LMSimple(obs,
                        model,
                        use_logpars=self['use_logpars'],
                        prior=self.search_prior,
                        lm_pars=self['lm_pars'])

        fitter.run_lm(guess)

        res=fitter.get_result()
        if res['flags']==0:
            print("        replacing cov")
            # try to replace the covariance
            h=1.0e-3
            m=5.0
            fitter.calc_cov(h, m)
            if res['flags'] != 0:
                print("        replacement failed")
                res['flags']=0

        return fitter


    def print_res(self,res):
        """
        print some stats
        """

        if 'neff' in res:
            mess="    s2n_w: %.1f  neff: %.1f  efficiency: %.2f"
            mess=mess % (res['s2n_w'],res['neff'],res['efficiency'])
        elif 'nfev' in res:
            mess="    s2n_w: %.1f  ntry: %d  nfev: %d"
            mess = mess % (res['s2n_w'],res['ntry'],res['nfev'])
        print(mess)

        print_pars(res['pars'],      front='        pars: ')
        print_pars(res['pars_err'],  front='        perr: ')
        print_pars(res['pars_true'], front='        true: ')


    def copy_to_output(self, res, i):
        super(NGMixSimISample,self).copy_to_output(res, i)
        d=self.data

        d['neff'][i] = res['neff']
        d['efficiency'][i] = res['efficiency']

    def get_dtype(self):
        dt=super(NGMixSimISample,self).get_dtype()
        dt += [('neff','f4'),
               ('efficiency','f4')]
        return dt


class NGMixSimISampleP(NGMixSimISample):
    def __init__(self, *args, **kw):
        super(NGMixSimISampleP,self).__init__(*args, **kw)

    def set_config(self, sim_conf, run_conf):
        from .util import get_shear_grid
        super(NGMixSimISampleP,self).set_config(sim_conf, run_conf)

        self.s1grid,self.s2grid = get_shear_grid(self.simc, self.conf)

    def _add_mcmc_stats(self, sampler):
        """
        measure prob of shear on a grid

        implementing equations 17-18 from B&A 2014
        """
        super(NGMixSimISampleP,self)._add_mcmc_stats(sampler)

        iweights = sampler.get_iweights()
        samples = sampler.get_samples()

        g1=samples[:,2]
        g2=samples[:,3]

        g_prior=self.g_prior

        if self.g_prior_during:
            prior_vals = g_prior.get_prob_array2d(g1,g2)
            w,=numpy.where(prior_vals > 0)
            if w.size == 0:
                raise TryAgainError("no prior vals > 0")
            iwsum = iweights[w].sum()
        else:
            iwsum = iweights.sum()

        lnp = numpy.zeros(self.conf['shear_grid']['dims']) - 9.999e9

        for i1,s1 in enumerate(self.s1grid):
            for i2,s2 in enumerate(self.s2grid):

                pjvals = g_prior.get_pj(g1, g2, s1, s2)

                if self.g_prior_during:
                    pjvals = pjvals[w]/prior_vals[w]

                p = (iweights[w]*pjvals).sum()/iwsum
                if p > 0:
                    lnp[i1,i2] = numpy.log(p)
                
        res=sampler.get_result()
        res['lnp_shear'] = lnp

    def copy_to_output(self, res, i):
        super(NGMixSimISampleP,self).copy_to_output(res, i)

        d=self.data
        d['lnp_shear'][i,:,:] = res['lnp_shear']

    def get_dtype(self):
        dt=super(NGMixSimISampleP,self).get_dtype()
        dims = self.conf['shear_grid']['dims']
        dt += [('lnp_shear','f8',dims)]
        return dt



class NGMixSimISampleComposite(NGMixSimISample):
    def run_max_fitter(self, imdict, model):
        """
        get the covariance matrix based on a max like fit

        raises if not found
        """
        from ngmix.fitting import FracdevFitter
        
        mess="model must be 'c' for composite, got '%s'" % model
        assert model=='composite',mess

        ipars=self['isample_pars']

        exp_fitter = super(NGMixSimISampleComposite,self).run_max_fitter(imdict, 'exp')
        dev_fitter = super(NGMixSimISampleComposite,self).run_max_fitter(imdict, 'dev')

        fres=self.run_fracdev(imdict, exp_fitter, dev_fitter)
        fracdev = fres['fracdev']

        TdByTe = self._get_TdByTe(exp_fitter, dev_fitter)

        for itry in xrange(1,ipars['max_tries']+1):
            max_guess=self.get_guess(imdict, guess_type='draw_truth')

            max_guess[4:] = log(max_guess[4:])

            fitter=self.run_composite(imdict, fracdev, TdByTe, max_guess)

            res=fitter.get_result()
            if res['flags']==0:
                break

        pars=res['pars']
        perr=res.get('pars_err',None)
        pars_cov=res.get('pars_cov',None)

        mess="    tries: %d  nfev: %d"
        mess=mess % (itry,res['nfev'])
        print(mess)
        print_pars(pars, front='    max pars: ')

        if res['flags'] != 0:
            if perr is None:
                print("    no cov")
            else:
                print("    did not converge")
            raise TryAgainError("bad max fit")

        res['fracdev'] = fres['fracdev']
        res['fracdev_err'] = fres['fracdev_err']

        print_pars(perr, front='    max perr: ')

        self.maxlike_res=res

        return fitter



    def _get_TdByTe(self, exp_fitter, dev_fitter):
        epars=exp_fitter.get_result()['pars']
        dpars=dev_fitter.get_result()['pars']

        if self['use_logpars']:
            Te = exp(epars[4])
            Td = exp(dpars[4])
        else:
            Te = epars[4]
            Td = dpars[4]
        TdByTe = Td/Te
        return TdByTe

    def run_fracdev(self, imdict, exp_fitter, dev_fitter):
        from ngmix.fitting import FracdevFitter
        ipars=self['isample_pars']

        epars=exp_fitter.get_result()['pars']
        dpars=dev_fitter.get_result()['pars']

        obs=imdict['obs']
        ffitter = FracdevFitter(obs, epars, dpars,
                                use_logpars=self['use_logpars'],
                                method=ipars['max_fitter'])
        for i in xrange(ipars['fracdev_tries']):
            ffitter.go(0.5 + 0.1*srandu())
            res=ffitter.get_result()
            if res['flags']==0:
                break

        if res['flags'] != 0:
            raise TryAgainError("failed fracdev fit")

        mess='    fracdev: %(fracdev).3f +/- %(fracdev_err).3f'
        mess = mess % res
        print(mess)

        return res

    def run_composite(self, imdict, fracdev, TdByTe, guess):
        """
        Fit the model to the galaxy

        we send some keywords so behavior can be different if this is a guess
        getter

        """
        from ngmix.fitting import LMComposite

        obs=imdict['obs']
        fitter=LMComposite(obs,
                           fracdev,
                           TdByTe,
                           use_logpars=self['use_logpars'],
                           prior=self.search_prior,
                           lm_pars=self['lm_pars'])

        fitter.go(guess)

        res=fitter.get_result()
        if res['flags']==0:
            print("        replacing cov")
            # try to replace the covariance
            h=1.0e-3
            m=5.0
            fitter.calc_cov(h, m)
            if res['flags'] != 0:
                print("        replacement failed")
                res['flags']=0

        return fitter

    def get_dtype(self):
        dt=super(NGMixSimISampleComposite,self).get_dtype()
        dt += [('fracdev','f4'),
               ('fracdev_err','f4')]
        return dt

    def copy_to_output(self, res, i):
        super(NGMixSimISampleComposite,self).copy_to_output(res, i)
        d=self.data

        mres=res['maxlike_res']
        d['fracdev'][i]     = mres['fracdev']
        d['fracdev_err'][i] = mres['fracdev_err']



def srandu(num=None):
    """
    Generate random numbers in the symmetric distribution [-1,1]
    """
    return 2*(numpy.random.random(num)-0.5)

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
    cmd='mv %s %s' % (local_file, output_file)
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

def check_g(g, maxval=0.97):
    gtot=sqrt(g[0]**2 + g[1]**2)
    if gtot > maxval:
        print("bad g:",gtot)
        raise TryAgainError("bad g")

def make_sheared_pars(pars, shear_g1, shear_g2):
    from ngmix import Shape
    shpars=pars.copy()

    sh=Shape(shpars[2], shpars[3])
    sh.shear(shear_g1, shear_g2)

    shpars[2]=sh.g1
    shpars[3]=sh.g2

    return shpars

def make_metacal_obs(psf_gmix, pars, model, noise_image, obs, nsub):
    """
    """
    from ngmix import Observation

    gm0=ngmix.GMixModel(pars, model)
    gm=gm0.convolve(psf_gmix)

    im = gm.make_image(obs.image.shape,
                       jacobian=obs.jacobian,
                       nsub=nsub)

    im += noise_image

    obs_new=Observation(im, jacobian=obs.jacobian, weight=obs.weight, psf=obs.psf)

    return obs_new


