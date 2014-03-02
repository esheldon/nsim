from __future__ import print_function

from sys import stderr
import time

import numpy
from numpy import array, zeros, cos, sin
from numpy.random import random as randu

def test(err_sigma, ngal, **keys):
    ts=ToySim(err_sigma, **keys)
    ts.go(ngal)
    return ts.get_data()

class ToySim(object):
    def __init__(self,
                 err_sigma,

                 ba_sigma=0.3,
                 shear=[0.08,0.0],

                 sampler='mcmc', # can also be 'true'

                 # for true sampler
                 nsamples=1000,

                 # for mcmc
                 nwalkers=80,
                 burnin=400,
                 nstep=200,
                 mca_a=2.0,
                 random_state=None,
                 print_step=1000,
                 make_plots=False):

        self.shear=array(shear)

        self.ba_sigma=ba_sigma
        self.err_sigma=err_sigma

        self.sampler_type=sampler

        self.nwalkers=nwalkers
        self.burnin=burnin
        self.nstep=nstep
        self.mca_a=mca_a

        self.nsamples=nsamples

        self.print_step=print_step
        self.make_plots=make_plots

        self.random_state=random_state

        self.set_distributions()

        if self.sampler_type=='mcmc':
            self._set_mcmc_sampler()

    def go(self, ngal):
        """
        Get P,Q,R for sampled points
        """

        self.make_output(ngal)

        i=0
        t0=time.time()
        for i in xrange(ngal):
            i1=i+1
            if (i1 % self.print_step)==0:
                print("gal %d/%d" % (i1, ngal), file=stderr)

            glike = self._get_galaxy_glike()

            P,Q,R,arate,neff=self._do_trials(glike)
            self._copy_to_output(i,P,Q,R,arate,neff)

        tm=time.time()-t0
        print("time per galaxy:",tm/ngal, file=stderr)
        print("total time:    ",tm, file=stderr)

    def get_data(self):
        """
        Get the results
        """
        if not hasattr(self,'data'):
            raise RuntimeError("you forgot to run go() !!")
        return self.data

    def _do_trials(self, glike):
        """
        Run the mcmc for a single object
        """
        if self.sampler_type=='mcmc':
            raise ValueError("make work for mcmc")
            runner=self._do_mcmc
        elif self.sampler_type=='true':
            runner=self._get_trials_true
        else:
            raise ValueError("bad sampler type: '%s'" % self.sampler_type)

        Pmax=0.0
        ntot=0

        Psum=0.0
        Qsum=numpy.zeros(2)
        Rsum=numpy.zeros( (2,2) )

        arate_sum=0.0

        niter=0
        neff=0.0
        while neff < self.nsamples:
            trials, arate = runner(glike)
            Pi, Qi, Ri = self._get_PQR_arrays(trials)

            Pmaxi = Pi.max()
            if Pmaxi > Pmax:
                Pmax=Pmaxi

            Psum += Pi.sum()
            Qsum += Qi.sum(axis=0)
            Rsum += Ri.sum(axis=0)


            neff = Psum/Pmax
            #print("    neff:",neff,file=stderr)

            arate_sum += arate
            niter+=1
            ntot += Pi.size

            if self.make_plots:
                self.do_make_plots(trials, glike)

        arate = arate_sum/niter

        P = Psum/ntot
        Q = Qsum/ntot
        R = Rsum/ntot
        return P,Q,R,arate,neff

    def _get_trials_true(self, glike):
        """
        Sample values directly from the actual like func
        """

        trials = zeros( (self.nsamples, 2) )
        trials[:,0], trials[:,1] = glike.sample(self.nsamples)

        arate=1.0
        return trials, arate


    def _get_PQR_arrays(self, trials):
        g1=trials[:,0]
        g2=trials[:,1]

        # expand around zero
        #Pi,Qi,Ri = self.g_prior.get_pqr(g1,g2)

        sh=self.shear
        Pi,Qi,Ri = self.g_prior.get_pqr_expand(g1,g2,sh[0],sh[1])

        return Pi, Qi, Ri

    def _get_PQR(self, trials):
        """
        get the marginalized P,Q,R from Bernstein & Armstrong

        If the prior is already in our mcmc chain, so we need to divide by the
        prior everywhere.

        zero prior values should be removed prior to calling
        """
        
        Pi, Qi, Ri = self._get_PQR_arrays(trials)

        neff=Pi.sum()/Pi.max()
        #Pi,Qi,Ri = self.g_prior.get_pqr_num(g1,g2)

        P = Pi.mean()
        Q = Qi.mean(axis=0)
        R = Ri.mean(axis=0)

        return P,Q,R,neff

    def _get_galaxy_glike(self):
        """
        Get the likelihood function
        """
        import ngmix

        g1,g2 = self.g_prior.sample2d(1)
        #g1,g2 = self.g_prior.sample2d_brute(1)

        shape=ngmix.shape.Shape(g1[0], g2[0])
        shape.shear( self.shear[0], self.shear[1] )

        #glike = self._get_glike(shape.g1, shape.g2)
        glike=self._add_error(shape)
        
        return glike

    def _add_error(self, shape):

        # first draw from gaussian centered at the true position
        glike0 = self._get_glike(shape.g1, shape.g2)

        g1p, g2p = glike0.sample()

        # now our likelihood is the same width but centered on the sampled
        # point
        return self._get_glike(g1p, g2p)

    def _get_glike(self, g1, g2):
        import ngmix
        glike = ngmix.priors.TruncatedGauss2D(g1,
                                              g2,
                                              self.err_sigma,
                                              self.err_sigma,
                                              1.0)

        return glike



    def _do_mcmc(self, glike):
        """
        Sample using mcmc
        """

        # set this global for the mcmc
        self.glike = glike

        sampler=self.mcmc_sampler

        sampler.reset()

        guess=self.get_mcmc_guess(glike)
        pos, prob, state = sampler.run_mcmc(guess, self.burnin)

        sampler.reset()
        pos, prob, state = sampler.run_mcmc(pos, self.nstep)

        trials  = sampler.flatchain
        arates = sampler.acceptance_fraction
        arate = arates.mean()

        return trials, arate


    def get_mcmc_guess(self, glike):
        """
        Draw guesses from the true likelihood
        """
        guess=zeros( (self.nwalkers,2) )
        guess[:,0],guess[:,1]=self.glike.sample(self.nwalkers)

        return guess

    def get_lnprob_mcmc(self, pars):
        """
        For the mcmc sampler
        """
        import ngmix
        from ngmix.gexceptions import GMixRangeError
        from ngmix.priors import LOWVAL

        g1,g2=pars[0],pars[1]

        lnp = self.glike.get_lnprob_nothrow(g1,g2)

        return lnp

    def _set_mcmc_sampler(self):
        """
        Instantiate the sampler
        """
        import emcee
        npars=2
        sampler = emcee.EnsembleSampler(self.nwalkers, 
                                        npars, 
                                        self.get_lnprob_mcmc,
                                        a=self.mca_a)
        if self.random_state is not None:

            # this is a property, runs set_state internally. sadly this will
            # fail silently which is the stupidest thing I have ever seen in my
            # entire life.  If I want to set the state it is important to me!
            
            print('    replacing random state',file=stderr)

            # OK, we will just hope that _random doesn't change names in the future.
            # but at least we get control back
            sampler._random = self.random_state

        self.mcmc_sampler = sampler

    def set_distributions(self):
        import ngmix
        self.g_prior=ngmix.priors.GPriorBA(self.ba_sigma)

    def do_make_plots(self, trials, glike):
        import mcmc
        tab=mcmc.plot_results(trials, names=['g1','g2'])

        if self.sampler_type != 'true':
            import biggles
            import esutil as eu
            # add some points drawn from true distribution
            sh=trials.shape
            true_trials=numpy.zeros( sh )
            true_trials[:,0],true_trials[:,1]=glike.sample(sh[0])

            

        key=raw_input('hit a key (q to quit): ')
        if key=='q':
            stop

    def _copy_to_output(self, i, P,Q,R, arate, neff):
        data=self.data
        data['P'][i] = P
        data['Q'][i,:] = Q
        data['R'][i,:,:] = R
        data['arate'][i] = arate
        data['neff'][i] = neff


    def make_output(self, ngal):
        dtype=[('P','f8'),
               ('Q','f8',2),
               ('R','f8', (2,2) ),
               ('arate','f8'),
               ('neff','f8')]

        self.data=zeros(ngal, dtype=dtype)


def quick(shear, err_sigma, ngal, ntrials):
    import ngmix
    P=numpy.zeros(ngal)
    Q=numpy.zeros( (ngal,2) )
    R=numpy.zeros( (ngal,2,2) )

    g_prior=ngmix.priors.GPriorBA(0.3)

    for i in xrange(ngal):
        i1=i+1
        if (i1 % 1000) == 0:
            print("%s/%s" % (i1,ngal),file=stderr)

        g1,g2=g_prior.sample2d(1)
        shape=ngmix.shape.Shape(g1[0],g2[0])
        shape.shear( shear[0], shear[1] )

        errdist=ngmix.priors.TruncatedGauss2D(shape.g1,shape.g2,
                                              err_sigma, err_sigma,
                                              1.0)
        g1p,g2p=errdist.sample()
        glike=ngmix.priors.TruncatedGauss2D(g1p,g2p,
                                            err_sigma, err_sigma,
                                            1.0)

        g1t,g2t=glike.sample(ntrials)
        Pi,Qi,Ri=g_prior.get_pqr_num(g1t,g2t,s1=shear[0],s2=shear[1])

        P[i]     = Pi.mean()
        Q[i,:]   = Qi.mean(axis=0)
        R[i,:,:] = Ri.mean(axis=0)

    return P,Q,R
