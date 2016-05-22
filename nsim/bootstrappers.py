from __future__ import print_function
from ngmix.bootstrap import Bootstrapper, MaxRunner
from ngmix.gexceptions import GMixRangeError, BootPSFFailure, BootGalFailure
from ngmix.guessers import TFluxGuesser, TFluxAndPriorGuesser

class MetacalMomentBootstrapper(Bootstrapper):

    def fit_max(self,
                gal_model,
                pars,
                prior=None,
                ntry=1,
                **kw):

        obs = self.mb_obs_list

        # note in this version we do not require  psfs

        guesser=self._get_max_guesser(prior=prior, **kw)

        runner=MaxRunner(obs, gal_model, pars, guesser,
                         prior=prior,
                         intpars=self.intpars,
                         use_logpars=self.use_logpars)

        runner.go(ntry=ntry)

        fitter=runner.fitter

        res=fitter.get_result()

        if res['flags'] != 0:
            raise BootGalFailure("failed to fit galaxy with maxlike")

        self.max_fitter=fitter



    def _get_max_guesser(self, prior=None, **kw):
        """
        this one we don't require psf was fit
        """

        Tguess=kw.get('Tguess',4.0)

        wsum=0.0
        fsum=0.0
        for obslist in self.mb_obs_list:
            for obs in obslist:
                wsum += obs.weight.sum()
                fsum += (obs.image*obs.weight).sum()
        if wsum > 0.0:
            flux_guess = fsum/wsum
        else:
            raise BootGalFailure("can't make a flux guess, weight is zero")

        if self.use_logpars:
            scaling='log'
        else:
            scaling='linear'

        if prior is None:
            guesser=TFluxGuesser(Tguess,
                                 flux_guess,
                                 scaling=scaling)
        else:
            guesser=TFluxAndPriorGuesser(Tguess,
                                         flux_guess,
                                         prior,
                                         scaling=scaling)
        return guesser




