from __future__ import print_function
import numpy
import ngmix
from ngmix.bootstrap import \
        Bootstrapper, \
        MaxRunner, \
        EMRunner, \
        MaxMetacalBootstrapper

from ngmix.gexceptions import GMixRangeError, BootPSFFailure, BootGalFailure
from ngmix.guessers import TFluxGuesser, TFluxAndPriorGuesser

class MetacalMomentBootstrapper(Bootstrapper):

    def get_em_fitter(self):
        return self.em_fitter

    def fit_em(self, Tguess, em_pars, ntry=1, **kw):

        obs = self.mb_obs_list

        runner=EMRunner(obs[0][0], Tguess, 1, em_pars)
        # note in this version we do not require  psfs

        runner.go(ntry=ntry)

        fitter=runner.fitter

        res=fitter.get_result()

        if res['flags'] != 0:
            raise BootGalFailure("failed to fit galaxy with maxlike")

        self.em_fitter=fitter



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




class MetaMomBootstrapper(Bootstrapper):
    def get_metamom_fitter(self):
        """
        get the maxlike fitter for the galaxy
        """
        if not hasattr(self,'metamom_fitter'):
            raise RuntimeError("you need to run fit_metamom successfully first")
        return self.metamom_fitter

    def get_fitter(self, fitter_type):
        """
        get fitter by name

        'max'
        'isample'
        etc.
        """
        if fitter_type=='metamom':
            fitter = self.get_metamom_fitter()
        else:
            fitter = super(MetaMomBootstrapper,self).get_fitter(fitter_type)
        return fitter

    def fit_metamom(self,
                    wt_gmix,
                    gal_model,
                    pars,
                    guess=None,
                    guess_widths=None,
                    prior=None,
                    extra_priors=None,
                    ntry=1):

        """
        fit the galaxy.  You must run fit_psf() successfully first
        """

        obs = self.mb_obs_list

        if not hasattr(self,'psf_flux_res'):
            self.fit_gal_psf_flux()

        guesser=self._get_max_guesser(
            guess=guess,
            prior=prior,
            widths=guess_widths,
        )

        runner=MetaMomRunner(
            wt_gmix,
            obs, gal_model, pars, guesser,
            prior=prior,
            intpars=self.intpars,
            use_logpars=self.use_logpars,
        )

        runner.go(ntry=ntry)

        fitter=runner.fitter

        res=fitter.get_result()

        if res['flags'] != 0:
            raise BootGalFailure("failed to fit galaxy with maxlike")

        self.metamom_fitter=fitter


class MetacalMetaMomBootstrapper(MaxMetacalBootstrapper):
    def fit_metacal(self,
                    wt_gmix,
                    psf_model,
                    gal_model,
                    pars,
                    psf_Tguess,
                    psf_fit_pars=None,
                    metacal_pars=None,
                    prior=None,
                    psf_ntry=5,
                    ntry=1,
                    **kw):
        """
        run metacalibration

        parameters
        ----------
        psf_model: string
            model to fit for psf
        gal_model: string
            model to fit
        pars: dict
            parameters for the maximum likelihood fitter
        psf_Tguess: float
            T guess for psf
        psf_fit_pars: dict
            parameters for psf fit
        metacal_pars: dict, optional
            Parameters for metacal, default {'step':0.01}
        prior: prior on parameters, optional
            Optional prior to apply
        psf_ntry: int, optional
            Number of times to retry psf fitting, default 5
        ntry: int, optional
            Number of times to retry fitting, default 1
        **kw:
            extra keywords for get_all_metacal
        """

        obs_dict = self._get_all_metacal(
            metacal_pars,
            **kw
        )

        res = self._do_metacal_max_fits(wt_gmix,
                                        obs_dict,
                                        psf_model, gal_model,
                                        pars, psf_Tguess,
                                        prior, psf_ntry, ntry,
                                        psf_fit_pars)

        self.metacal_res = res

    def _do_metacal_max_fits(self, wt_gmix,obs_dict, psf_model, gal_model, pars, 
                             psf_Tguess, prior, psf_ntry, ntry, 
                             psf_fit_pars):

        # overall flags, or'ed from each bootstrapper
        res={'mcal_flags':0}
        for key in sorted(obs_dict):
            # run a regular Bootstrapper on these observations
            boot = MetaMomBootstrapper(
                obs_dict[key],
                use_logpars=self.use_logpars,
                intpars=self.intpars,
                use_round_T=self.use_round_T,
                find_cen=self.find_cen,
                verbose=self.verbose,
            )

            boot.fit_psfs(psf_model, psf_Tguess, ntry=psf_ntry,
                          fit_pars=psf_fit_pars,
                          skip_already_done=False)

            boot.fit_metamom(
                wt_gmix,
                gal_model,
                pars,
                prior=prior,
                ntry=ntry,
            )

            boot.set_round_s2n(fitter_type='metamom')

            tres=boot.get_metamom_fitter().get_result()
            rres=boot.get_round_result()

            res['mcal_flags'] |= tres['flags']
            res['mcal_flags'] |= rres['flags']

            tres['s2n_r'] = rres['s2n_r']
            tres['T_r'] = rres['T_r']
            tres['psf_T_r'] = rres['psf_T_r']

            gpsf_sum = numpy.zeros(2)
            Tpsf_sum = 0.0
            npsf=0
            for obslist in boot.mb_obs_list:
                for obs in obslist:
                    if hasattr(obs,'psf_nopix'):
                        #print("    summing nopix")
                        g1,g2,T=obs.psf_nopix.gmix.get_g1g2T()
                    else:
                        g1,g2,T=obs.psf.gmix.get_g1g2T()
                    gpsf_sum[0] += g1
                    gpsf_sum[1] += g2
                    Tpsf_sum += T
                    npsf+=1

            tres['gpsf'] = gpsf_sum/npsf
            tres['Tpsf'] = Tpsf_sum/npsf

            res[key] = tres

        return res

class MetaMomRunner(MaxRunner):
    """
    wrapper to generate guesses and run the fitter a few times
    """
    def __init__(self, wt_gmix, *args, **kw):
        super(MetaMomRunner,self).__init__(*args, **kw)
        self.wt_gmix=wt_gmix

    def _go_lm(self, ntry=1):

        if self.intpars is not None:
            npoints=self.intpars['npoints']
        else:
            npoints=None

        fitclass=self._get_lm_fitter_class()

        for i in xrange(ntry):
            guess=self.guesser()
            fitter=fitclass(self.obs,
                            self.model,
                            self.wt_gmix,
                            lm_pars=self.send_pars,
                            use_logpars=self.use_logpars,
                            use_round_T=self.use_round_T,
                            npoints=npoints,
                            prior=self.prior)

            fitter.go(guess)

            res=fitter.get_result()
            if res['flags']==0:
                break

        res['ntry'] = i+1
        self.fitter=fitter

    def _get_lm_fitter_class(self):
        from .fitting import LMSimple
        return LMSimple



    def _get_lm_fitter_class(self):
        return ngmix.fitting.LMMetaMomSimple


