"""
Simulate images and fit them, currently MCMC only

Currently importing ngmix locally because it is so slow; numba does not yet
cache object files. Will make global import when caching is part of numba

additional dependence on
    emcee for MCMC fitting and
    fitsio if checkpointing

example sim config file in yaml format
--------------------------------------
name: "nsim-dg01"

psf_model: "gauss"
psf_T: 4.0
psf_shape: [0.0, 0.0]

obj_model: "dev"
obj_T_mean: 16.0
obj_T_sigma_frac: 0.3

obj_counts_mean: 100.0
obj_counts_sigma_frac: 0.3

shear: [0.01,0.0]

nsub: 16

example run config file in yaml format
--------------------------------------
run: "ngmix-dg01r33"
sim: "nsim-dg01"

fit_model: "dev"

nwalkers: 40
burnin:   400
nstep:    200
mca_a:    3.0

# we normalize splits by split for is2n==0
desired_err: 2.0e-05
nsplit0: 60000

s2n_vals: [ 15, 21, 30, 42, 60, 86, 122, 174, 247, 352, 500] 

"""

import os
from sys import stderr
import time
import pprint

import numpy
from numpy.random import random as randu
from numpy.random import randn


# region over which to render images and calculate likelihoods
NSIGMA_RENDER=5.0

# minutes
#DEFAULT_CHECKPOINTS=[5,30,60,100]
DEFAULT_CHECKPOINTS=[30,90]
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
        import ngmix

        # seed a new MT random generator from devrand
        # and return the object
        self.random_state=get_random_state_devrand()

        self.set_config(sim_conf, run_conf)
        self.update(self.conf)
        self.update(keys)

        self.shear=self.simc['shear']
        self.nsub=self.simc['nsub']
        self.nsigma_render=self.simc.get('nsigma_render',NSIGMA_RENDER)

        self.check_pqr_shear()

        self.s2n=s2n
        self.npairs=npairs

        self.ring=self.get('ring',True)

        self.obj_model=self.simc['obj_model']
        self.fit_model=self['fit_model']
        self.npars=ngmix.gmix.get_model_npars(self.fit_model)

        self.make_plots=keys.get('make_plots',False)
        self.plot_base=keys.get('plot_base',None)

        self.setup_checkpoints(**keys)

        if self.data is None:
            self.make_struct()

        self.set_priors()
        self.make_psf()
        self.set_noise()

        pprint.pprint(self, stream=stderr)

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

        if self['expand_shear_true']:
            self.shear_expand=self.shear
            print >>stderr,'nsim: expanding about shear:',self.shear_expand
        else:
            self.shear_expand=None

    def get_data(self):
        """
        Get a ref to the data array with the fit results
        """
        return self.data

    def run_sim(self):
        """
        Run the simulation, fitting psf and all pairs
        """
        self.fit_psf()

        self.start_timer()

        i=0
        npairs=self.npairs
        for ipair in xrange(npairs):
            print >>stderr,'%s/%s' % (ipair+1,npairs)

            self.ipair=ipair
            if self.data['processed'][i]:
                i += 2 # skip the pair
            else:
                while True:
                    try:
                        reslist=self.process_pair()
                        break
                    except TryAgainError as err:
                        print >>stderr,str(err)

                self.copy_to_output(reslist[0], i)
                i += 1
                self.copy_to_output(reslist[1], i)
                i += 1

            self.set_elapsed_time()
            self.try_checkpoint()

        self.set_elapsed_time()
        print >>stderr,'time minutes:',self.tm_minutes
        print >>stderr,'time per pair sec:',self.tm/npairs
        print >>stderr,'time per image sec:',self.tm/(2*npairs)

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
        import ngmix

        imdicts = self.get_noisy_image_pair()
        reslist=[]
        for key in imdicts:

            fitter=self.fit_galaxy(imdicts[key])
            res=fitter.get_result()
            if res['flags'] != 0:
                raise TryAgainError("failed at %s" % key)

            reslist.append(res)
            self.print_res(res,imdicts[key]['pars'])

            if self.make_plots:
                self.do_make_plots(fitter,key)

        return reslist

    def do_make_plots(self, fitter, key):
        """
        Write a plot file of the trials
        """
        import biggles
        biggles.configure('default','fontsize_min',1.0)

        if 'isample' in self['fitter']:
            dims=(800,1100)
            width,height=800,1100
        else:
            width,height=1100,1100

        tp=fitter.make_plots(title=self.fit_model)
        if isinstance(tp, tuple):
            p,wp=tp

            trials_plot=self.plot_base+'-%06d-%s-trials.png' % (self.ipair,key)
            trials_wplot=self.plot_base+'-%06d-%s-wtrials.png' % (self.ipair,key)

            print >>stderr,trials_plot
            p.write_img(width,height,trials_plot)
            print >>stderr,trials_wplot
            wp.write_img(width,height,trials_wplot)

        else:
            trials_plot=self.plot_base+'-%06d-%s-trials.png' % (self.ipair,key)
            print >>stderr,trials_plot
            tp.write_img(width,height,trials_plot)


    def fit_galaxy(self, imdict):
        """
        Fit according to the requested method
        """

        fitter_type=self['fitter']
        #if fitter_type == 'mcmc':
        if 'mcmc' in fitter_type:
            fitter = self.fit_galaxy_mcmc(imdict)
        elif fitter_type == 'lm':
            fitter = self.fit_galaxy_lm(imdict)

        elif fitter_type == 'mh':
            fitter = self.fit_galaxy_mh(imdict)

        elif fitter_type == 'isample':
            fitter = self.fit_galaxy_isample(imdict)

        elif fitter_type == 'isample-adapt':
            fitter = self.fit_galaxy_isample_adapt(imdict)

        elif fitter_type=='isample-anze':
            fitter = self.fit_galaxy_isample_anze(imdict)

        else:
            raise ValueError("bad fitter type: '%s'" % fitter_type)

        return fitter

    def fit_galaxy_isample(self, imdict):
        """
        Fit the model to the galaxy using important sampling
        """
        import ngmix

        # do checking for guess here?
        if self['guess_type']=='truth':
            sampler=self.get_simple_proposal_truth(imdict)
            trials,ln_probs=sampler.sample(self['n_samples'])
        else:
            raise ValueError("support guess type: '%s'" % self['guess_type'])

        fitter=ngmix.fitting.ISampleSimple(imdict['image'],
                                           imdict['wt'],
                                           imdict['jacobian'],
                                           self.fit_model,
                                           trials,
                                           ln_probs,


                                           # not optional
                                           cen_prior=self.cen_prior,
                                           g_prior=self.g_prior,
                                           T_prior=self.T_prior,
                                           counts_prior=self.counts_prior,

                                           shear_expand=self.shear_expand,

                                           psf=self.psf_gmix_fit,

                                           do_pqr=True,
                                           do_lensfit=True)
        fitter.go()
        res = fitter.get_result()

        ew=res['eff_iweight']
        print >>stderr,'    eff iweight:',ew,'eff samples:',res['eff_n_samples']

        if False:
            pprint.pprint(res)
        return fitter


    def fit_galaxy_isample_adapt(self, imdict):
        """
        Fit the model to the galaxy using important sampling
        """
        import ngmix

        # this loop only makes sense for random guesses
        # do checking for guess here?
        if self['proposal_type']=='truth':
            # note this suffers from bad "sigma" values
            sampler=self.get_simple_proposal_truth(imdict)
        elif self['proposal_type']=='maxlike':
            # suffers from potential bad covariance matrices
            sampler=self.get_simple_proposal_maxlike(imdict)
        elif self['proposal_type']=='mcmc':
            # suffers from potential bad covariance matrices
            sampler=self.get_simple_proposal_mcmc(imdict)
        else:
            raise ValueError("support guess type: '%s'" % self['guess_type'])

        fitter=ngmix.fitting.ISampleSimpleAdapt(imdict['image'],
                                                imdict['wt'],
                                                imdict['jacobian'],
                                                self.fit_model,
                                                sampler, # starting sampler component

                                                n_per=self['n_per'],
                                                max_components=self['max_components'],
                                                max_fdiff=self['max_fdiff'],

                                                # not optional
                                                cen_prior=self.cen_prior,
                                                g_prior=self.g_prior,
                                                T_prior=self.T_prior,
                                                counts_prior=self.counts_prior,

                                                g_prior_during=self['g_prior_during'],

                                                shear_expand=self.shear_expand,

                                                psf=self.psf_gmix_fit,

                                                do_pqr=True,
                                                do_lensfit=True)
        fitter.go()
        res = fitter.get_result()

        if not fitter.has_converged():
            raise TryAgainError("not converged")

        return fitter



    def get_simple_proposal_maxlike(self, imdict):
        """
        These guesses are better for low s/n but still pretty bad
        """
        import ngmix

        if self['guess_type'] not in ['draw_priors','truth_random']:
            raise ValueError("use draw_priors or truth_random for lm start")

        ntry=10
        for i in xrange(ntry):
            fitter=self.fit_galaxy_lm(imdict)
            res=fitter.get_result()
            if res['flags']==0:
                break

        if res['flags'] != 0:
            raise TryAgainError("failed fit LM after %s tries" % ntry)

        pars=res['pars']
        perr=res['pars_err']

        # soften covariance, need to tune
        #perr += numpy.array([0.005,0.005,0.005,0.005,0.1,1.0])

        cen_dist    = ngmix.priors.Student2D(pars[0], pars[1], perr[0], perr[1])
        g_dist      = ngmix.priors.TruncatedStudentPolar(pars[2], pars[3], perr[2], perr[3], 1.0)
        #T_dist      = ngmix.priors.LogNormal(pars[4], perr[4])
        #counts_dist = ngmix.priors.LogNormal(pars[5], perr[5])
        T_dist      = ngmix.priors.StudentPositive(pars[4], perr[4])
        counts_dist = ngmix.priors.StudentPositive(pars[5], perr[5])
 
        sampler=SimpleSampler(cen_dist, g_dist, T_dist, counts_dist)
        return sampler

    def get_simple_proposal_mcmc(self, imdict):
        """
        These guesses are better for low s/n but still pretty bad
        """
        import ngmix

        fitter=self.fit_galaxy_mcmc(imdict)

        res=fitter.get_result()
        pars=res['pars']
        perr=res['pars_err']*1.2

        cen_dist    = ngmix.priors.Student2D(pars[0], pars[1], perr[0], perr[1])
        g_dist      = ngmix.priors.TruncatedStudentPolar(pars[2], pars[3], perr[2], perr[3], 1.0)
        T_dist      = ngmix.priors.LogNormal(pars[4], perr[4])
        counts_dist = ngmix.priors.LogNormal(pars[5], perr[5])
 
        sampler=SimpleSampler(cen_dist, g_dist, T_dist, counts_dist)
        return sampler


    def get_simple_proposal_truth(self, imdict):
        """
        These guesses are better for low s/n but still pretty bad
        """
        import ngmix

        pars=imdict['pars']

        cen0=pars[0]
        cen1=pars[1]
        # should tune based on a max-like fit or something
        g1=pars[2]
        g2=pars[3]
        T = pars[4]
        counts = pars[5]

        cen_sigma=self.simc['cen_sigma']
        g_sigma=0.2
        T_sigma = self.simc['obj_T_sigma_frac']*T
        counts_sigma = self.simc['obj_counts_sigma_frac']*counts

        #cen_dist=ngmix.priors.CenPrior(cen0, cen1, cen_sigma, cen_sigma)
        #g_dist=ngmix.priors.TruncatedGaussianPolar(g1, g2, g_sigma, g_sigma, 1.0)
        #T_dist=ngmix.priors.LogNormal(T, T_sigma)
        #counts_dist=ngmix.priors.LogNormal(counts, counts_sigma)

        cen_dist    = ngmix.priors.Student2D(cen0, cen1, cen_sigma, cen_sigma)
        g_dist      = ngmix.priors.TruncatedStudentPolar(g1, g2, g_sigma, g_sigma, 1.0)
        T_dist=ngmix.priors.LogNormal(T, T_sigma)
        counts_dist=ngmix.priors.LogNormal(counts, counts_sigma)
 
        sampler=SimpleSampler(cen_dist, g_dist, T_dist, counts_dist)
        return sampler


    def fit_galaxy_isample_anze(self, imdict):
        """
        Fit the model to the galaxy using anze's game sampler
        """
        import ngmix

        
        # faking this for now
        nwalkers=1

        full_guess=imdict['pars'].copy()

        fitter=ngmix.fitting.ISampleBDFAnze(imdict['image'],
                                            imdict['wt'],
                                            imdict['jacobian'],

                                            cen_prior=self.cen_prior,
                                            g_prior=self.g_prior,
                                            T_prior=self.T_prior,
                                            counts_prior=self.counts_prior,

                                            bfrac_prior=self.bfrac_prior,

                                            g_prior_during=self['g_prior_during'],

                                            full_guess=full_guess,

                                            shear_expand=self.shear_expand,

                                            psf=self.psf_gmix_fit,
                                            nwalkers=nwalkers,
                                            nstep=self['nstep'],
                                            burnin=self['burnin'],
                                            do_pqr=True,
                                            do_lensfit=True)

        fitter.go()
        res = fitter.get_result()
        if False:
            pprint.pprint(res)

        return res




    def fit_galaxy_isample_old(self, imdict):
        """
        Fit the model to the galaxy using important sampling
        """
        import ngmix

        trials,ln_probs=self.draw_isamples_priors()
        #trials,ln_probs=self.draw_isamples_from_true(imdict['pars'])

        fitter=ngmix.fitting.ISampleSimple(imdict['image'],
                                           imdict['wt'],
                                           imdict['jacobian'],
                                           self.fit_model,

                                           trials,
                                           ln_probs,

                                           cen_prior=self.cen_prior,
                                           g_prior=self.g_prior,
                                           T_prior=self.T_prior,
                                           counts_prior=self.counts_prior,

                                           shear_expand=self.shear_expand,

                                           psf=self.psf_gmix_fit,

                                           n_samples=self['n_samples'],
                                           do_pqr=True,
                                           do_lensfit=True)
        fitter.go()
        res = fitter.get_result()
        if False:
            import pprint
            pprint.pprint(res)

        ew=res['eff_iweight']
        print >>stderr,'    eff iweight:',ew,'eff samples:',res['eff_n_samples']

        return res

    def draw_isamples_from_true(self, pars):
        """
        Draw samples for importance sampling

        Forcing certain types of functions here
        """

        import ngmix

        if self.npars != 6:
            raise ValueError("support guess from non-simple!")

        n_samples=self['n_samples']
        trials = numpy.zeros( (n_samples, self.npars) )
        ln_probs = numpy.zeros(n_samples)

        cen_sigma=self.simc['cen_sigma']
        g_sigma=0.3
        T_sigma = self.simc['obj_T_sigma_frac']*pars[4]
        counts_sigma = self.simc['obj_counts_sigma_frac']*pars[5]

        cen_dist=ngmix.priors.CenPrior(pars[0], pars[1], cen_sigma, cen_sigma)
        g_dist = ngmix.priors.CenPrior(pars[2],pars[3], g_sigma, g_sigma)
        T_dist=ngmix.priors.LogNormal(pars[4], T_sigma)
        counts_dist=ngmix.priors.LogNormal(pars[5], counts_sigma)


        trials[:,0],trials[:,1] = cen_dist.sample(n=n_samples)
        trials[:,2],trials[:,3] = g_dist.sample(n=n_samples)
        trials[:,4]             = T_dist.sample(nrand=n_samples)
        trials[:,5]             = counts_dist.sample(nrand=n_samples)

        for i in xrange(n_samples):
            pars=trials[i,:]

            lnp = 0.0

            lnp += cen_dist.get_lnprob(pars[0], pars[1])
            lnp += g_dist.get_lnprob(pars[2], pars[3])
            lnp += T_dist.get_lnprob_scalar(pars[4])
            lnp += counts_dist.get_lnprob_scalar(pars[5])

            ln_probs[i] = lnp

        return trials, ln_probs


    def draw_isamples_priors(self):
        """
        Draw samples for importance sampling

        For now just draw from priors, but we will need to change this
        """

        import ngmix

        if self.npars != 6:
            raise ValueError("support isample guess from non-simple!")

        n_samples=self['n_samples']
        trials = numpy.zeros( (n_samples, self.npars) )
        ln_probs = numpy.zeros(n_samples)

        trials[:,0],trials[:,1]=self.cen_prior.sample(n=n_samples)
        trials[:,2],trials[:,3]=self.g_prior.sample2d(n_samples)
        trials[:,4]=self.T_prior.sample(nrand=n_samples)
        trials[:,5]=self.counts_prior.sample(nrand=n_samples)

        cen_prior=self.cen_prior
        g_prior=self.g_prior
        T_prior=self.T_prior
        counts_prior=self.counts_prior
        for i in xrange(n_samples):
            pars=trials[i,:]

            lnp = 0.0

            lnp += cen_prior.get_lnprob(pars[0], pars[1])
            lnp += g_prior.get_lnprob_scalar2d(pars[2], pars[3])
            lnp += self.T_prior.get_lnprob_scalar(pars[4])
            lnp += self.counts_prior.get_lnprob_scalar(pars[5])

            ln_probs[i] = lnp

        return trials, ln_probs

    def fit_galaxy_mcmc(self, imdict):
        """
        Fit the model to the galaxy
        """
        import ngmix

        if self.fit_model=='bdf':
            if self['restart']:
                fitter=self.run_bdf_mcmc_fitter_with_restart(imdict)
            else:
                fitter=self.run_bdf_mcmc_fitter(imdict)
        else:
            fitter=self.run_simple_mcmc_fitter(imdict)

        return fitter


    def fit_galaxy_mh(self, imdict):
        """
        Fit the model to the galaxy
        """
        if self.fit_model in ['gauss','exp','dev']:
            fitter=self.run_simple_mh_fitter(imdict)
        else:
            raise ValueError("implement other mh fitters")
        return fitter


    def get_guess(self, imdict, n=1):
        """
        Get a full guess, nwalkers x npars
        """
        guess_type=self['guess_type']
        if guess_type=='draw_truth':
            print >>stderr,'guessing randomized truth'
            full_guess=self.get_guess_from_pars(imdict['pars'], n=n)
        elif guess_type=='draw_priors':
            full_guess=self.get_guess_draw_priors(n=n)
        elif guess_type=='draw_maxlike':
            full_guess=self.get_guess_draw_maxlike(imdict, n=n)
        else:
            raise ValueError("bad guess type: '%s'" % guess_type)
        
        return full_guess

    def get_mcmc_fitter(self, imdict, full_guess):
        """
        Get the appropriate MCMC fitter
        """
        import ngmix

        return fitter

    def run_bdf_mcmc_fitter_with_restart(self, imdict):
        """
        Get a bdf (Bulge-Disk Fixed size ratio) mcmc fitter
        """
        import ngmix
        from ngmix.fitting import LOW_ARATE

        nwalkers = self['nwalkers']
        ntry=self['ntry']

        for ipass in [0,1]:
            if ipass==0:
                print >>stderr,'    pass 1'
                full_guess=self.get_guess(imdict, n=nwalkers)
            else:
                best_pars=fitter.best_pars
                ngmix.fitting.print_pars(best_pars,
                                         front='    pass 2 pars: ',
                                         stream=stderr)
                full_guess=self.get_guess_from_pars(best_pars, n=nwalkers)

            min_arate=self['min_arate'][ipass]
            for itry in xrange(ntry):
                fitter=ngmix.fitting.MCMCBDF(imdict['image'],
                                             imdict['wt'],
                                             imdict['jacobian'],

                                             cen_prior=self.cen_prior,
                                             g_prior=self.g_prior,
                                             T_prior=self.T_prior,
                                             counts_prior=self.counts_prior,

                                             bfrac_prior=self.bfrac_prior,

                                             g_prior_during=self['g_prior_during'],

                                             full_guess=full_guess,

                                             shear_expand=self.shear_expand,

                                             psf=self.psf_gmix_fit,
                                             nwalkers=nwalkers,
                                             nstep=self['nstep'],
                                             burnin=self['burnin'],
                                             mca_a=self['mca_a'][ipass],
                                             ntry=1,
                                             min_arate=min_arate,
                                             random_state=self.random_state,
                                             do_pqr=True,
                                             do_lensfit=True)
                fitter.go()
                if fitter.arate < min_arate and itry < (ntry-1):
                    nwalkers = nwalkers*2
                    print >>stderr,'        trying nwalkers:',nwalkers
                    full_guess=self.get_guess(imdict, n=nwalkers)
                else:
                    break


        return fitter

    def run_bdf_mcmc_fitter(self, imdict):
        """
        Get a bdf (Bulge-Disk Fixed size ratio) mcmc fitter
        """
        import ngmix
        from ngmix.fitting import LOW_ARATE

        nwalkers = self['nwalkers']

        full_guess=self.get_guess(imdict, n=nwalkers)

        ntry=self['ntry']
        min_arate=self['min_arate']

        for i in xrange(ntry):
            fitter=ngmix.fitting.MCMCBDF(imdict['image'],
                                         imdict['wt'],
                                         imdict['jacobian'],

                                         cen_prior=self.cen_prior,
                                         g_prior=self.g_prior,
                                         T_prior=self.T_prior,
                                         counts_prior=self.counts_prior,

                                         bfrac_prior=self.bfrac_prior,

                                         g_prior_during=self['g_prior_during'],

                                         full_guess=full_guess,

                                         shear_expand=self.shear_expand,

                                         psf=self.psf_gmix_fit,
                                         nwalkers=nwalkers,
                                         nstep=self['nstep'],
                                         burnin=self['burnin'],
                                         mca_a=self['mca_a'],
                                         #ntry=ntry,
                                         min_arate=min_arate,
                                         random_state=self.random_state,
                                         do_pqr=True,
                                         do_lensfit=True)
            fitter.go()
            res=fitter.get_result()

            arate=res['arate']
            if arate >= min_arate:
                break
            elif i < (ntry-1):
                print >>stderr,"        arate",arate,"<",min_arate
                nwalkers_old=nwalkers
                nwalkers = nwalkers*2
                print >>stderr,'        trying nwalkers:',nwalkers

                #full_guess=self.get_guess(imdict, n=nwalkers)

                # start where we left off
                #sampler=fitter.get_sampler()
                w=fitter.lnprobs.argmax()
                best_pars=fitter.trials[w, :]
                full_guess=self.get_guess_from_pars(best_pars, n=nwalkers)

                #full_guess=sampler.chain[:, -2:, :].reshape(nwalkers,self.npars)

                #full_guess=numpy.zeros( (nwalkers,self.npars) )
                #chain=sampler.chain
                #pos=fitter.get_last_pos()
                #full_guess[0:nwalkers_old]=pos
                #full_guess[nwalkers_old:]=pos
                #full_guess[0:nwalkers_old]=chain[:, -1:, :]
                #full_guess[nwalkers_old:]=chain[:, -2:, :]
                #full_guess[:,:] = trials[-nwalkers:, :]

                print '       ',full_guess.shape

        if arate < min_arate:
            if self['keep_low_arate']:
                fitter._result['flags'] =0
            else:
                fitter._result['flags'] |= LOW_ARATE


        # trying keeping low arate ones
        return fitter

    def run_simple_mcmc_fitter(self, imdict):
        """
        Get a bdf (Bulge-Disk Fixed size ratio) mcmc fitter
        """
        import ngmix

        full_guess=self.get_guess(imdict, n=self['nwalkers'])

        fitter=ngmix.fitting.MCMCSimple(imdict['image'],
                                        imdict['wt'],
                                        imdict['jacobian'],
                                        self.fit_model,

                                        cen_prior=self.cen_prior,
                                        g_prior=self.g_prior,
                                        T_prior=self.T_prior,
                                        counts_prior=self.counts_prior,

                                        g_prior_during=self['g_prior_during'],

                                        full_guess=full_guess,

                                        shear_expand=self.shear_expand,

                                        psf=self.psf_gmix_fit,
                                        nwalkers=self['nwalkers'],
                                        nstep=self['nstep'],
                                        burnin=self['burnin'],
                                        mca_a=self['mca_a'],
                                        random_state=self.random_state,
                                        do_pqr=True,
                                        do_lensfit=True)
        fitter.go()
        return fitter

    def run_simple_mh_fitter(self, imdict):
        """
        Run metropolis hastings
        """
        import ngmix

        # get guess and step size from lm
        ntry=10
        for i in xrange(ntry):
            lm_fitter=self.fit_galaxy_lm(imdict)
            res=lm_fitter.get_result()
            if res['flags']==0:
                break

        if res['flags'] != 0:
            raise TryAgainError("failed fit LM after %s tries" % ntry)

        guess=res['pars']
        step_sizes = 0.5*res['pars_err'].copy()
        # soften step sizes, need to tune.
        #step_sizes += numpy.array([0.001, 0.001, 0.005, 0.005, 0.1, 0.1])

        fitter=ngmix.fitting.MHSimple(imdict['image'],
                                      imdict['wt'],
                                      imdict['jacobian'],
                                      self.fit_model,

                                      cen_prior=self.cen_prior,
                                      g_prior=self.g_prior,
                                      T_prior=self.T_prior,
                                      counts_prior=self.counts_prior,

                                      g_prior_during=self['g_prior_during'],

                                      shear_expand=self.shear_expand,

                                      psf=self.psf_gmix_fit,

                                      nstep=self['nstep'],
                                      burnin=self['burnin'],
                                      guess=guess,
                                      step_sizes=step_sizes,
                                      min_arate=self['min_arate'],
                                      ntry=self['ntry'],

                                      do_pqr=True,
                                      do_lensfit=True)
        fitter.go()
        return fitter


    def get_guess_draw_priors(self, n=1):
        """
        Get a guess drawn from the priors

        assume simple for now
        """
        if self.fit_model == 'bdf':
            guess=self.get_guess_draw_priors_bdf(n)
        else:
            guess=self.get_guess_draw_priors_simple(n)

        if n==1:
            guess=guess[0,:]

        return guess


    def get_guess_draw_priors_simple(self, n=1):
        """
        Get a guess drawn from the priors
        """
        import ngmix

        #print >>stderr,'    ** guessing from priors'
        if self.npars != 6:
            raise ValueError("wrong npars for simple!")

        guess=numpy.zeros( (n, self.npars) )

        if self.cen_prior is not None:
            guess[:,0],guess[:,1]=self.cen_prior.sample(n=n)
        guess[:,2],guess[:,3]=self.g_prior.sample2d(n)
        guess[:,4]=self.T_prior.sample(nrand=n)
        guess[:,5]=self.counts_prior.sample(nrand=n)

        return guess

    def get_guess_draw_priors_bdf(self, n=1):
        """
        Get a guess drawn from the priors
        """
        import ngmix

        #print >>stderr,'    guessing from priors'
        if self.npars != 7:
            raise ValueError("wrong npars for bdf!")

        guess=numpy.zeros( (n, self.npars) )

        guess[:,0],guess[:,1]=self.cen_prior.sample(n=n)
        guess[:,2],guess[:,3]=self.g_prior.sample2d(n)
        guess[:,4]=self.T_prior.sample(nrand=n)
        
        # sample from the bfrac distribution
        counts=self.counts_prior.sample(nrand=n)
        #bfracs = self.bfrac_prior.sample(n)

        # uniform
        # make sure there are some at high and low bfrac
        bfracs = randu(n)
        if n > 2:
            bfracs=numpy.zeros(n)
            nhalf=n/2
            bfracs[0:nhalf] = 0.01*randu(nhalf)
            bfracs[nhalf:] = 0.99+0.01*randu(nhalf)
        else:
            bfracs = randu(n)
        dfracs = 1.0-bfracs

        # bulge flux
        guess[:,5] = bfracs*counts
        # disk flux
        guess[:,6] = dfracs*counts

        return guess

    def get_guess_draw_maxlike(self, imdict, n=1):
        """
        Get the maximum likelihood fit and draw from that using
        width from the fit

        start lm near the truth
        """
        import ngmix

        ntry=10
        for i in xrange(ntry):
            lm_guess=self.get_guess_from_pars(imdict['pars'], n=1)
            fitter=self.fit_galaxy_lm(imdict, guess=lm_guess)
            res=fitter.get_result()
            if res['flags']==0:
                break

        if res['flags'] != 0:
            raise TryAgainError("failed fit LM after %s tries" % ntry)

        pars=res['pars']
        perr=res['pars_err']
        ngmix.fitting.print_pars(pars, front='        lmpars: ',stream=stderr)
        ngmix.fitting.print_pars(perr, front='        lmperr: ',stream=stderr)

        guess=self.get_guess_from_pars(pars, n=n, width=perr)

        return guess

    def get_guess_from_pars(self, pars, n=1, width=None):
        """
        Get a guess centered on the truth
        """
        if self.fit_model=='bdf':
            guess=self.get_guess_from_pars_bdf(pars,n=n,width=width)
        else:
            guess=self.get_guess_from_pars_simple(pars,n=n,width=width)

        if n==1:
            guess=guess[0,:]

        return guess

    def get_guess_from_pars_bdf(self, pars, n=1, width=None):
        """
        Get a guess centered on the truth

        width is relative for T and counts
        """

        print >>stderr,'guessing pars bdf'
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
        guess[:,4] = self.get_positive_guess(pars[4],n,width=width[4])
        guess[:,5] = self.get_positive_guess(pars[5],n,width=width[5])
        guess[:,6] = self.get_positive_guess(pars[6],n,width=width[6])

        return guess


    def get_guess_from_pars_simple(self, pars, n=1, width=None):
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
        guess[:,4] = self.get_positive_guess(pars[4],n,width=width[4])
        guess[:,5] = self.get_positive_guess(pars[5],n,width=width[5])

        return guess

    def get_shape_guess(self, g1, g2, n, width=[0.01,0.01]):
        """
        Get guess, making sure in range
        """
        import ngmix
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

    def get_positive_guess(self, val, n, width=0.01):
        """
        Get guess, making sure positive
        """
        from ngmix.gexceptions import GMixRangeError

        if val <= 0.0:
            raise GMixRangeError("val <= 0: %s" % val)

        vals=numpy.zeros(n)-9999.0
        while True:
            w,=numpy.where(vals <= 0)
            if w.size == 0:
                break
            else:
                vals[w] = val*(1.0 + width*srandu(w.size))

        return vals


    def fit_galaxy_lm(self, imdict, guess=None):
        """
        Fit the model to the galaxy
        """
        import ngmix

        if guess is None:
            guess=self.get_lm_guess(imdict)

        im=imdict['image']
        wt=imdict['wt']
        j=imdict['jacobian']
        psf=self.psf_gmix_fit
        
        counts_prior = self.counts_prior

        fitter=ngmix.fitting.LMSimple(im, wt, j,
                                      self.fit_model,
                                      guess,

                                      lm_pars=self['lm_pars'],

                                      cen_prior=self.cen_prior,
                                      g_prior=self.g_prior,
                                      T_prior=self.T_prior,
                                      counts_prior=counts_prior,

                                      psf=psf)

        fitter.go()

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

    def print_res(self,res,true_pars):
        """
        print some stats
        """
        import ngmix

        if 'arate' in res:
            print >>stderr,'    arate:',res['arate'],'s2n_w:',\
                    res['s2n_w'],'nuse:',res['nuse']
        elif 'nfev' in res:
            print >>stderr,'    nfev:',res['nfev'],'s2n_w:',res['s2n_w']

        ngmix.fitting.print_pars(true_pars, front='    true: ',stream=stderr)
        ngmix.fitting.print_pars(res['pars'],front='    pars: ',stream=stderr)
        ngmix.fitting.print_pars(res['pars_err'],front='    perr: ',stream=stderr)

    def fit_psf(self):
        """
        Fit the pixelized psf to a model

        PSF

        """
        import ngmix
        from ngmix.gexceptions import GMixMaxIterEM

        print >>stderr,'fitting psf'
        imsky,sky=ngmix.em.prep_image(self.psf_image)

        em=ngmix.em.GMixEM(imsky)

        tol=self.get('psf_tol',1.0e-5)
        maxiter=self.get('psf_maxiter',1000)

        while True:
            guess=self.get_psf_guess()
            print >>stderr,'psf guess:'
            print >>stderr,guess
            try:
                em.go(guess, sky, tol=tol,maxiter=maxiter)
                break
            except GMixMaxIterEM as e:
                print >>stderr,str(e)
                print >>stderr,'re-trying'

        self.psf_gmix_fit=em.get_gmix()
        print >>stderr,'psf fit:'
        print >>stderr,self.psf_gmix_fit

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

    def set_priors(self):
        """
        Set all the priors
        """
        import ngmix

        print >>stderr,"setting priors"
        T=self.simc['obj_T_mean']
        T_sigma = self.simc['obj_T_sigma_frac']*T
        counts=self.simc['obj_counts_mean']
        counts_sigma = self.simc['obj_counts_sigma_frac']*counts

        cen_sigma=self.simc['cen_sigma']
        self.cen_prior=ngmix.priors.CenPrior(0.0, 0.0, cen_sigma, cen_sigma)

        self.g_prior=ngmix.priors.GPriorBA(0.3)
        self.T_prior=ngmix.priors.LogNormal(T, T_sigma)
        self.counts_prior=ngmix.priors.LogNormal(counts, counts_sigma)

        self.bfrac_prior=ngmix.priors.BFrac()

    def make_psf(self):
        """
        make the psf gaussian mixture model
        """
        import ngmix

        print >>stderr,"making psf"

        self.psf_dims, self.psf_cen=self.get_dims_cen(self.simc['psf_T'])

        pars=[self.psf_cen[0],
              self.psf_cen[1],
              self.simc['psf_shape'][0],
              self.simc['psf_shape'][1],
              self.simc['psf_T'],
              1.0]
        self.psf_gmix_true=ngmix.gmix.GMixModel(pars, self.simc['psf_model'])
        
        self.psf_image=self.psf_gmix_true.make_image(self.psf_dims,
                                                     nsub=self.nsub)
    
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
        

        print >>stderr,"setting noise"

        imdict=self.get_image_pair(random=False)
        im=imdict['im1']['image']
        skysig2 = (im**2).sum()/self.s2n**2
        skysig = numpy.sqrt(skysig2)

        s2n_check = numpy.sqrt( (im**2).sum()/skysig**2 )
        print >>stderr,"S/N goal:",self.s2n,"found:",s2n_check

        self.skysig=skysig
        self.ivar=1.0/skysig**2


    def get_noisy_image_pair(self):
        """
        Get an image pair, with noise added
        """
        imdict=self.get_image_pair()
        self.add_noise(imdict['im1']['image'])
        self.add_noise(imdict['im2']['image'])

        wt=numpy.zeros(imdict['im1']['image'].shape) + self.ivar
        imdict['im1']['wt']=wt
        imdict['im2']['wt']=wt
        return imdict

    def add_noise(self, im):
        """
        Add gaussian random noise
        """

        im[:,:] += self.skysig*randn(im.size).reshape(im.shape)

    def get_image_pair(self, random=True):
        """
        get a model image

        If random is True, use draw random values from the priors.
        Otherwise use the mean of the priors

        """
        import ngmix

        pars1, pars2, cen_offset = self.get_pair_pars(random=random)

        gm1_pre=ngmix.gmix.GMixModel(pars1, self.obj_model)
        gm2_pre=ngmix.gmix.GMixModel(pars2, self.obj_model)

        gm1  = gm1_pre.convolve(self.psf_gmix_true)
        gm2  = gm2_pre.convolve(self.psf_gmix_true)

        T = gm1.get_T()
        dims, cen = self.get_dims_cen(T)

        # jacobian is at center before offset so the prior
        # will center on "zero"
        j=ngmix.jacobian.UnitJacobian(cen[0], cen[1])

        cen[0] += cen_offset[0]
        cen[1] += cen_offset[1]

        gm1.set_cen(cen[0], cen[1])
        gm2.set_cen(cen[0], cen[1])

        nsub = self.nsub
        im1=gm1.make_image(dims, nsub=nsub)
        im2=gm2.make_image(dims, nsub=nsub)

        pars_true1=numpy.array(pars1)
        pars_true2=numpy.array(pars2)
        pars_true1[0] += cen_offset[0]
        pars_true1[1] += cen_offset[1]
        pars_true2[0] += cen_offset[0]
        pars_true2[1] += cen_offset[1]

        out={'im1':{'pars':pars_true1,'gm_pre':gm1_pre,'gm':gm1,'image':im1,'jacobian':j},
             'im2':{'pars':pars_true2,'gm_pre':gm2_pre,'gm':gm2,'image':im2,'jacobian':j}}
        return out

    def get_pair_pars(self, random=True):
        """
        Get pair parameters

        """
        import ngmix


        if random:
            cen_offset=self.cen_prior.sample()

            g = self.g_prior.sample1d(1)
            g=g[0]
            rangle1 = randu()*2*numpy.pi
            if self.ring:
                rangle2 = rangle1 + numpy.pi/2.0
            else:
                rangle2 = randu()*2*numpy.pi
            g1_1 = g*numpy.cos(2*rangle1)
            g2_1 = g*numpy.sin(2*rangle1)
            g1_2 = g*numpy.cos(2*rangle2)
            g2_2 = g*numpy.sin(2*rangle2)

            T=self.T_prior.sample()
            counts=self.counts_prior.sample()
        else:
            cen_offset=[0.0, 0.0]
            g1_1=0.0
            g2_1=0.0
            g1_2=0.0
            g2_2=0.0
            T=self.T_prior.mean
            counts=self.counts_prior.mean

        shape1=ngmix.shape.Shape(g1_1, g2_1)
        shape2=ngmix.shape.Shape(g1_2, g2_2)

        shear=self.shear
        shape1.shear(shear[0], shear[1])
        shape2.shear(shear[0], shear[1])

        # center is just placeholder for now
        pars1=[0.0, 0.0, shape1.g1, shape1.g2, T]
        pars2=[0.0, 0.0, shape2.g1, shape2.g2, T]

        if self.obj_model=='bdf':
            bfrac=self.bfrac_prior.sample()
            counts_b=bfrac*counts
            counts_d=(1.0-bfrac)*counts
            pars1+=[counts_b, counts_d]
            pars2+=[counts_b, counts_d]
        else:
            pars1+=[counts]
            pars2+=[counts]

        return pars1, pars2, cen_offset

    def get_dims_cen(self, T):
        """
        Based on T, get the required dimensions and a center
        """
        sigma=numpy.sqrt(T/2.)
        dims = [2.*sigma*self.nsigma_render]*2
        cen = [(dims[0]-1.)/2.]*2

        return dims, cen

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

        print >>stderr,'checkpointing at',self.tm_minutes,'minutes'
        success=write_fits(self.checkpoint_file, self.data)


    def copy_to_output(self, res, i):
        """
        Copy results into the output
        """
        d=self.data
        d['processed'][i] = 1
        d['pars'][i,:] = res['pars']
        d['pcov'][i,:,:] = res['pars_cov']

        if 'P' in res:
            d['P'][i] = res['P']
            d['Q'][i,:] = res['Q']
            d['R'][i,:,:] = res['R']
            d['g'][i,:] = res['g']
            d['gsens'][i,:] = res['g_sens']
            d['nuse'][i] = res['nuse']
        else:
            d['nfev'][i] = res['nfev']

    def make_struct(self):
        """
        Make the output array
        """
        npars=self.npars

        dt=[('processed','i2'),
            ('pars','f8',npars),
            ('pcov','f8',(npars,npars))]

        if 'lm' in self['fitter']:
            dt += [('nfev','i4')]
        else:
            dt += [('P','f8'),
                   ('Q','f8',2),
                   ('R','f8',(2,2)),
                   ('g','f8',2),
                   ('gsens','f8',2),
                   ('nuse','i4')]
        self.data=numpy.zeros(self.npairs*2, dtype=dt)

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
    print >>stderr,"writing local file:",local_file

    with fitsio.FITS(local_file,'rw',clobber=True) as fobj:
        fobj.write(data)

    if local_file==output_file:
        return

    # remove if it exists
    try:
        os.remove(output_file)
    except:
        pass

    # try a few times
    print >>stderr,"moving to:",output_file
    cmd='mv -v %s %s' % (local_file, output_file)
    for i in xrange(5):
        stat=os.system(cmd)
        if stat==0:
            print >>stderr,'success'
            success=True
            break
        else:
            print >>stderr,'error moving file, trying again'
            time.sleep(5)
            success=False

    return success

def get_random_state_devrand():
    """
    Seed the numpy random state from /dev/random
    """
    seed = get_devrand_uint()
    print >>stderr,'seed from devrand:',seed
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
