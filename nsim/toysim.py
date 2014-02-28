from __future__ import print_function

from sys import stderr
import time

import numpy
from numpy import array, zeros, cos, sin
from numpy.random import random as randu

def test(err_sigma, npair, **keys):
    ts=ToySim(err_sigma, **keys)
    ts.go(npair)
    return ts.get_data()

class ToySim(object):
    def __init__(self,
                 err_sigma,

                 ba_sigma=0.3,
                 shear=[0.01,0.0],

                 sampler='mcmc', # can also allow truth

                 # for true sampler
                 nsamples_true=2000,

                 # for mcmc
                 nwalkers=80,
                 burnin=400,
                 nstep=200,
                 mca_a=2.0,
                 random_state=None,
                 print_step=1,
                 make_plots=False):

        self.shear=array(shear)

        self.ba_sigma=ba_sigma
        self.err_sigma=err_sigma

        self.sampler_type=sampler

        self.nwalkers=nwalkers
        self.burnin=burnin
        self.nstep=nstep
        self.mca_a=mca_a

        self.nsamples_true=nsamples_true

        self.print_step=print_step
        self.make_plots=make_plots

        self.random_state=random_state

        self.set_distributions()

        if self.sampler_type=='mcmc':
            self._set_mcmc_sampler()

    def go(self, npair):
        """
        Run an MCMC for each trial
        """

        self.make_output(npair)

        i=0
        t0=time.time()
        for ipair in xrange(npair):
            ipair1=ipair+1
            if (ipair1 % self.print_step)==0:
                print("pair %d/%d" % (ipair1, npair), file=stderr)
            glikes_pair = self.get_pair()

            for glike in glikes_pair:
                P,Q,R,arate=self._do_trials(glike)
                self._copy_to_output(i,P,Q,R,arate)

                i+=1

        tm=time.time()-t0
        print("time per image:",tm/npair/2., file=stderr)
        print("time per pair: ",tm/npair, file=stderr)
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
            trials, arate = self._do_mcmc(glike)
        elif self.sampler_type=='true':
            trials = self._get_trials_true(glike)
            arate=1.0
        else:
            raise ValueError("bad sampler type: '%s'" % self.sampler_type)

        if self.make_plots:
            self.do_make_plots(trials)

        P,Q,R = self._get_PQR(trials)
        return P,Q,R,arate

    def _get_trials_true(self, glike):
        """
        Sample values directly from the actual like func
        """

        trials = zeros( (self.nsamples_true, 2) )
        trials[:,0], trials[:,1] = glike.sample(self.nsamples_true)

        return trials

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
        #try:
        #    lnp = self.glike.get_lnprob(g1,g2)
        #except GMixRangeError:
        #    lnp=LOWVAL

        return lnp

    def _get_PQR(self, trials):
        """
        get the marginalized P,Q,R from Bernstein & Armstrong

        If the prior is already in our mcmc chain, so we need to divide by the
        prior everywhere.

        zero prior values should be removed prior to calling
        """
        
        g1=trials[:,0]
        g2=trials[:,1]

        # expand around zero
        Pi,Qi,Ri = self.g_prior.get_pqr(g1,g2)

        P = Pi.mean()
        Q = Qi.mean(axis=0)
        R = Ri.mean(axis=0)

        return P,Q,R

    def _get_mean_pqr(self, Pi, Qi, Ri):
        """
        Get the mean P,Q,R marginalized over priors.  Optionally weighted for
        importance sampling
        """


        return P,Q,R


    def get_pair(self):
        """
        Get g1,g2 values with error added
        """

        g = self.g_prior.sample1d(1)
        g=g[0]

        rangle_first = randu()*2*numpy.pi
        rangle_sec = rangle_first + numpy.pi/2.0

        glike_first = self._get_glike_from_true(g, rangle_first)
        glike_sec   = self._get_glike_from_true(g, rangle_sec)

        return glike_first, glike_sec

    def _get_glike_from_true(self, g, theta):
        import ngmix
        g1 = g*cos(2*theta)
        g2 = g*sin(2*theta)

        shape=ngmix.shape.Shape(g1, g2)
        shape.shear( self.shear[0], self.shear[1] )

        g1_true=shape.g1
        g2_true=shape.g2

        glike=self._add_error(g1_true, g2_true)
        
        return glike

    def _add_error(self, g1, g2):

        # first draw from gaussian centered at the true position
        glike = self._get_glike(g1, g2)

        g1p, g2p = glike.sample()

        # now our likelihood is the same width but centered on the one with error
        new_glike = self._get_glike(g1p, g2p)

        return new_glike

    def _get_glike(self, g1, g2):
        import ngmix
        glike = ngmix.priors.TruncatedGauss2D(g1,
                                              g2,
                                              self.err_sigma,
                                              self.err_sigma,
                                              1.0)

        return glike

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

    def do_make_plots(self, trials):
        import mcmc
        mcmc.plot_results(trials, names=['g1','g2'])
        key=raw_input('hit a key (q to quit): ')
        if key=='q':
            stop

    def _copy_to_output(self, i, P,Q,R,arate):
        data=self.data
        data['P'][i] = P
        data['Q'][i,:] = Q
        data['R'][i,:,:] = R
        data['arate'][i] = arate


    def make_output(self, npair):
        dtype=[('P','f8'),
               ('Q','f8',2),
               ('R','f8', (2,2) ),
               ('arate','f8')]

        ntot = npair*2
        self.data=zeros(ntot, dtype=dtype)


