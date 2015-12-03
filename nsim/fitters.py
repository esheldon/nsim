"""
fit simulated images
"""
from __future__ import print_function
import os
import time
from pprint import pprint

import numpy
from numpy import array, zeros, ones, log, log10, exp, sqrt, diag
from numpy import where, isfinite
from numpy.random import random as randu
from numpy.random import randn

import ngmix
from ngmix import srandu
from ngmix.fitting import print_pars
from ngmix.gexceptions import GMixRangeError
from ngmix.observation import Observation
from ngmix.gexceptions import GMixMaxIterEM
from ngmix.shape import Shape
from ngmix.gmix import GMixModel

from ngmix.gexceptions import BootPSFFailure, BootGalFailure

from .sim import NGMixSim
from .util import write_fits, TryAgainError, load_gmixnd
from . import files

# minutes
DEFAULT_CHECKPOINTS=[30,60,90,110]

class FitterBase(dict):
    def __init__(self, sim, run_conf, ngal, **keys):

        self.sim=sim
        self._setup(run_conf, **keys)

        self['ngal']=ngal

        self._set_prior()

        self._setup_checkpoints(**keys)
        if self.data is None:
            self._make_struct()

        pprint(self)

    def go(self):
        """
        process the requested number of simulated galaxies
        """

        self._start_timer()
        self.tm_sim=0.0
        self.tm_fit=0.0

        nprocessed=0
        for igal in xrange(self['ngal']):
            print('%s/%s' % (igal+1,self['ngal']) )

            if self.data['processed'][igal]:
                continue

            while True:
                nprocessed += 1
                try:

                    tm0=time.time()
                    imdict = self.sim.get_image()
                    self.tm_sim += time.time()-tm0

                    tm0=time.time()
                    res=self.process_one(imdict)
                    self.tm_fit += time.time()-tm0

                    res['shear_true'] = imdict['gal_info']['shear']
                    res['shear_index'] = imdict['gal_info']['shear_index']
                    self._copy_to_output(res, igal)

                    self._set_elapsed_time()
                    self._try_checkpoint()

                    break
                except TryAgainError as err:
                    print(str(err))


        self._set_elapsed_time()

        print("nprocessed (including failures):",nprocessed)
        print('time minutes:',self.tm_minutes)
        print('time per image sec:',self.tm/nprocessed)
        print('time to simulate:',self.tm_sim/nprocessed)
        print('time to fit:',self.tm_fit/nprocessed)

    def process_one(self, imdict):
        """
        perform the fit
        """

        print("   ",imdict['obs'].image.shape)

        fitter=self._dofit(imdict)

        res=fitter.get_result()
        if res['flags'] != 0:
            raise TryAgainError("failed at %s "
                                "flags %s" % (key,res['flags']))

        res['pars_true'] = imdict['pars']
        res['model_true'] = imdict['model']
        res['s2n_true'] = imdict['s2n']

        self._print_res(res)

        if self['make_plots']:
            self._make_plots(fitter,key)

        return res

    def get_data(self):
        """
        Get a ref to the data array with the fit results
        """
        return self.data

    def _make_plots(self, fitter, key):
        """
        Write a plot file of the trials
        """
        raise RuntimeError("over-ride me")

    def _dofit(self, imdict):
        """
        Fit according to the requested method
        """

        raise RuntimeError("over-ride me")

    def _copy_to_output(self, res, i):
        d=self.data

        d['processed'][i] = 1
        d['model_true'][i] = res['model_true']
        d['s2n_true'][i] = res['s2n_true']
        d['pars_true'][i,:] = res['pars_true']

        d['shear_true'][i] = res['shear_true'].g1,res['shear_true'].g2
        d['shear_index'][i] = res['shear_index']

    def _setup(self, run_conf, **keys):
        """
        Check and set the configurations
        """

        self['nrand'] = self.get('nrand',1)
        if self['nrand'] is None:
            self['nrand']=1

        if self.sim['name'] != run_conf['sim']:
            err="sim name in run config '%s' doesn't match sim name '%s'"
            raise ValueError(err % (run_conf['sim'],self.sim['name']))

        self.update(run_conf)

        self['use_logpars']=self.get('use_logpars',False)
        self['npars']=ngmix.gmix.get_model_npars(self['fit_model'])

        if isinstance(self.sim['obj_model'], dict):
            npars=2
        else:
            npars = ngmix.gmix.get_model_npars(self.sim['obj_model'])
        self['npars_true']=npars

        self['make_plots']=keys.get('make_plots',False)
        self['plot_base']=keys.get('plot_base',None)

    def _set_prior(self):
        raise RuntimeError("over-ride me")

    def _start_timer(self):
        """
        Set the elapsed time so far
        """
        self.tm0 = time.time()

    def _set_elapsed_time(self):
        """
        Set the elapsed time so far
        """

        self.tm = time.time()-self.tm0
        self.tm_minutes = self.tm/60.0


    def _setup_checkpoints(self, **keys):
        """
        Set up checkpoint times, file, and sent data
        """

        self.checkpoints     = keys.get('checkpoints',DEFAULT_CHECKPOINTS)
        self.n_checkpoint    = len(self.checkpoints)
        self.checkpointed    = [False]*self.n_checkpoint

        self.checkpoint_file=keys.get('checkpoint_file',None)
        self._set_checkpoint_data(**keys)

        if self.checkpoint_file is not None:
            self.do_checkpoint=True
        else:
            self.do_checkpoint=False


    def _set_checkpoint_data(self, **keys):
        """
        Look for checkpoint data, file etc.
        """
        self.data=None

        checkpoint_data=keys.get('checkpoint_data',None)
        if checkpoint_data is not None:
            self.data=checkpoint_data


    def _try_checkpoint(self):
        """
        If we should make a checkpoint, do so
        """

        should_checkpoint, icheck = self._should_checkpoint()

        if should_checkpoint:
            self._write_checkpoint()
            self.checkpointed[icheck]=True


    def _should_checkpoint(self):
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


    def _write_checkpoint(self):
        """
        Write the checkpoint file

        The file is written to the cwd and if this is not the final
        destination, an attempt is made to move it there.  This may fail and if
        so a message is printed.

        """

        print('checkpointing at',self.tm_minutes,'minutes')
        success=write_fits(self.checkpoint_file, self.data)


    def _get_dtype(self):
        """
        Set generic dtype
        """

        dt=[('processed','i2'),
            ('model_true','S3'),
            ('s2n_true','f8'),
            ('pars_true','f8',self['npars_true']),
            ('shear_true','f8',2),
            ('shear_index','i4'),
           ]

        return dt

    def _make_struct(self):
        """
        Make the output array
        """

        dt=self._get_dtype()
        self.data=numpy.zeros(self['ngal'], dtype=dt)

class SimpleFitterBase(FitterBase):

    def _copy_to_output(self, res, i):

        super(SimpleFitterBase,self)._copy_to_output(res, i)

        d=self.data

        d['pars'][i,:] = res['pars']
        d['pars_cov'][i,:,:] = res['pars_cov']

        d['g'][i,:] = res['g']
        d['g_cov'][i,:,:] = res['g_cov']

        d['s2n_w'][i] = res['s2n_w']

    def _set_prior(self):
        """
        Set all the priors
        """
        import ngmix
        from ngmix.joint_prior import PriorSimpleSep

        ppars=self['priors']

        gp = ppars['g']
        if gp['type']=='ba':
            self.g_prior = ngmix.priors.GPriorBA(gp['sigma'])
        elif gp['type']=='cosmos':
            self.g_prior=ngmix.priors.make_gprior_cosmos_sersic(type='erf')
        elif gp['type']=='great-des':
            self.g_prior = ngmix.priors.GPriorGreatDES(pars=gp['pars'],
                                                       gmax=1.0)

        else:
            raise ValueError("implement other g prior")

        self.g_prior_during = self.get('g_prior_during',False)
        if self.g_prior_during:
            print("    Using full g prior for fits")
            fit_g_prior = self.g_prior
        else:
            print("    Using flat g prior for fits")
            fit_g_prior = ngmix.priors.ZDisk2D(1.0)

        print("using input search prior")

        Tp = ppars['T']
        if Tp['type']=="truth":
            print("using true T pdf for prior")
            if self['use_logpars']:

                print("    converting to log")
                T=self.sim['obj_T_mean']
                T_sigma = self.sim['obj_T_sigma_frac']*T
                logT_mean, logT_sigma=ngmix.priors.lognorm_convert(T,T_sigma)
                fit_T_prior = ngmix.priors.Normal(logT_mean, logT_sigma)

            else:
                fit_T_prior = self.sim.T_pdf

        elif Tp['type']=="gmixnd":
            fit_T_prior = load_gmixnd(Tp)

        elif Tp['type']=='normal':
            Tpars=Tp['pars']
            fit_T_prior=ngmix.priors.Normal(Tpars[0], Tpars[1])

        elif Tp['type']=='lognormal':

            if self['use_logpars']:
                print("    converting T prior to log")
                logT_mean, logT_sigma=ngmix.priors.lognorm_convert(Tp['mean'],
                                                                   Tp['sigma'])
                fit_T_prior = ngmix.priors.Normal(logT_mean, logT_sigma)

            else:
                fit_T_prior = ngmix.priors.LogNormal(Tp['mean'],Tp['sigma'])

        elif Tp['type']=="two-sided-erf":
            T_prior_pars = Tp['pars']
            fit_T_prior=ngmix.priors.TwoSidedErf(*T_prior_pars)
        else:
            raise ValueError("bad Tprior: '%s'" % Tp['type'])

        cp=ppars['counts']
        if cp['type']=="truth":
            print("using true counts pdf for prior")

            if self['use_logpars']:

                print("    converting to log")
                if self.sim['simulator']=="galsim":
                    fluxspec=self.sim['obj_model']['flux']
                    counts       = fluxspec['mean']
                    counts_sigma = fluxspec['sigma']
                else:
                    counts       = self.sim['obj_counts_mean']
                    counts_sigma = self.sim['obj_counts_sigma_frac']*counts

                logc_mean, logc_sigma=ngmix.priors.lognorm_convert(counts,
                                                                   counts_sigma)
                fit_counts_prior = ngmix.priors.Normal(logc_mean, logc_sigma)

            else:
                fit_counts_prior = self.sim.counts_pdf

        elif cp['type']=="gmixnd":

            fit_counts_prior = load_gmixnd(cp)

        elif cp['type']=='lognormal':

            if self['use_logpars']:
                print("    converting counts prior to log")
                logcounts_mean, logcounts_sigma=ngmix.priors.lognorm_convert(cp['mean'],
                                                                   cp['sigma'])
                fit_counts_prior = ngmix.priors.Normal(logcounts_mean, logcounts_sigma)

            else:
                fit_counts_prior = ngmix.priors.LogNormal(cp['mean'],cp['sigma'])


        elif cp['type']=='normal':
            cpars=cp['pars']
            fit_counts_prior=ngmix.priors.Normal(cpars[0], cpars[1])

        elif cp['type']=="two-sided-erf":
            counts_prior_pars = cp['pars']
            fit_counts_prior=ngmix.priors.TwoSidedErf(*counts_prior_pars)
        else:
            raise ValueError("bad counts prior: '%s'" % cp['type'])

        cp=ppars['cen']
        if cp['type']=="truth":
            fit_cen_prior=self.sim.cen_pdf
        elif cp['type'] == "normal2d":
            fit_cen_sigma=cp['sigma']
            fit_cen_prior=ngmix.priors.CenPrior(0.0,
                                                0.0,
                                                fit_cen_sigma,
                                                fit_cen_sigma)
        else:
            raise ValueError("bad cen prior: '%s'" % cp['type'])


        self.prior = PriorSimpleSep(fit_cen_prior,
                                    fit_g_prior,
                                    fit_T_prior,
                                    fit_counts_prior)



    def _get_dtype(self):
        """
        get the dtype for the output struct
        """
        npars=self['npars']

        dt=super(SimpleFitterBase,self)._get_dtype()
        dt += [
            ('pars','f8',npars),
            ('pars_cov','f8',(npars,npars)),
            ('g','f8',2),
            ('g_cov','f8',(2,2)),
            ('s2n_w','f8')
        ]

        return dt

class MaxFitter(SimpleFitterBase):

    def _dofit(self, imdict):
        """
        Fit according to the requested method
        """

        boot=ngmix.Bootstrapper(imdict['obs'],
                                use_logpars=self['use_logpars'],
                                verbose=False)

        Tguess=self.sim.get('psf_T',4.0)
        ppars=self['psf_pars']
        try:
            boot.fit_psfs(ppars['model'], Tguess, ntry=ppars['ntry'])
        except BootPSFFailure:
            raise TryAgainError("failed to fit psf")

        mconf=self['max_pars']
        try:
            boot.fit_max(self['fit_model'],
                         mconf['pars'],
                         prior=self.prior,
                         ntry=mconf['ntry'])

            if mconf['pars']['method']=='lm':
                boot.try_replace_cov(mconf['cov_pars'])
        except BootGalFailure:
            raise TryAgainError("failed to fit galaxy")

        fitter=boot.get_max_fitter() 

        return fitter

    def _print_res(self,res):
        """
        print some stats
        """

        if 'nfev' in res:
            mess="    s2n_w: %.1f  ntry: %d  nfev: %d"
            mess = mess % (res['s2n_w'],res['ntry'],res['nfev'])
            print(mess)

        print_pars(res['pars'],      front='        pars: ')
        print_pars(res['pars_err'],  front='        perr: ')
        print_pars(res['pars_true'], front='        true: ')

    def _make_plots(self, fitter, key):
        """
        Write a plot file of the trials
        """
        import biggles

        biggles.configure('default','fontsize_min',1.0)

        width,height=1100,1100

        pdict=fitter.make_plots(do_residual=True, title=self.fit_model)

        resid_pname=self.plot_base+'-%06d-%s-gal-resid.png' % (self.ipair,key)
        print(resid_pname)
        rplt=fitter.plot_residuals()
        rplt[0][0].write_img(width,height,resid_pname)
    
        '''
        copy into ones that have samples
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
        '''


    def _get_dtype(self):
        """
        get the dtype for the output struct
        """
        npars=self['npars']

        dt=super(MaxFitter,self)._get_dtype()
        dt += [
            ('nfev','i4'),
            ('ntry','i4')
        ]
        return dt


    def _copy_to_output(self, res, i):

        super(MaxFitter,self)._copy_to_output(res, i)

        d=self.data

        if 'nfev' in res:
            d['nfev'][i] = res['nfev']
            # set outside of fitter
            d['ntry'][i] = res['ntry']


class MaxMetacalFitter(MaxFitter):
    """
    metacal with a maximum likelihood fit
    """
    def _dofit(self, imdict):
        """
        Fit according to the requested method
        """

        self.imdict=imdict


        obs=imdict['obs']
        mdict = self._do_fits(obs)
        res=mdict['res']
        fitter=mdict['fitter']

        self.fitter=fitter

        return fitter

    def _do_fits(self, obs):
        """
        the basic fitter for this class
        """
        from ngmix import Bootstrapper

        intpars=self.get('intpars',None) 

        boot=Bootstrapper(obs,
                          use_logpars=self['use_logpars'],
                          intpars=intpars,
                          verbose=False)

        Tguess=self.sim.get('psf_T',4.0)
        ppars=self['psf_pars']

        psf_fit_pars = ppars.get('fit_pars',None)
    
        try:
            boot.fit_psfs(ppars['model'], Tguess, ntry=ppars['ntry'],fit_pars=psf_fit_pars)
        except BootPSFFailure:
            raise TryAgainError("failed to fit psf")

        mconf=self['max_pars']

        try:
            boot.fit_max(self['fit_model'],
                         mconf['pars'],
                         prior=self.prior,
                         ntry=mconf['ntry'])


            if mconf['pars']['method']=='lm':
                boot.try_replace_cov(mconf['cov_pars'])

            self._do_metacal(boot)

        except BootPSFFailure:
            raise TryAgainError("failed to fit metacal psfs")
        except BootGalFailure:
            raise TryAgainError("failed to fit galaxy")

        fitter=boot.get_max_fitter() 
        res=fitter.get_result()

        mres = boot.get_metacal_max_result()
        res.update(mres)

        return {'fitter':fitter,
                'boot':boot,
                'res':res}

    def _do_metacal(self, boot, metacal_obs=None):

        ppars=self['psf_pars']
        mconf=self['max_pars']
        Tguess=self.sim.get('psf_T',4.0)
        psf_fit_pars = ppars.get('fit_pars',None)

        boot.fit_metacal_max(
            ppars['model'],
            self['fit_model'],
            mconf['pars'],
            Tguess,
            psf_fit_pars=psf_fit_pars,
            prior=self.prior,
            ntry=mconf['ntry'],
            metacal_pars=self['metacal_pars'],
            metacal_obs=metacal_obs
        )


    def _print_res(self,res):
        """
        print some stats
        """

        super(MaxMetacalFitter,self)._print_res(res)
        print("    mcal s2n_r:",res['mcal_s2n_r'],'mcal s2n_simple:',res['mcal_s2n_simple'])
        print_pars(res['mcal_pars'],       front='    mcal pars: ')
        print_pars(res['mcal_R'].ravel(),  front='    mcal R:    ')
        print_pars(res['mcal_Rpsf'],       front='    mcal Rpsf: ')

        print("    mcal c:",res['mcal_c'][0], res['mcal_c'][1])

    def _get_dtype(self):
        """
        get the dtype for the output struct
        """
        dt=super(MaxMetacalFitter,self)._get_dtype()

        npars=self['npars']
        dt += [
            ('mcal_pars','f8',npars),
            ('mcal_pars_cov','f8',(npars,npars)),
            ('mcal_g','f8',2),
            ('mcal_g_cov','f8', (2,2) ),
            ('mcal_g_noshear','f8',2),
            ('mcal_c','f8',2),
            ('mcal_s2n_r','f8'),
            ('mcal_s2n_simple','f8'),
            ('mcal_R','f8',(2,2)),
            ('mcal_Rpsf','f8',2),
            ('mcal_gpsf','f8',2),
        ]
        return dt


    def _copy_to_output(self, res, i):
        """
        copy parameters specific to this class
        """
        super(MaxMetacalFitter,self)._copy_to_output(res, i)

        d=self.data

        d['mcal_pars'][i] = res['mcal_pars']
        d['mcal_pars_cov'][i] = res['mcal_pars_cov']
        d['mcal_g'][i] = res['mcal_g']
        d['mcal_g_cov'][i] = res['mcal_g_cov']
        d['mcal_s2n_r'][i] = res['mcal_s2n_r']
        d['mcal_s2n_simple'][i] = res['mcal_s2n_simple']

        d['mcal_g_noshear'][i] = res['mcal_pars_noshear'][2:2+2]
        d['mcal_c'][i] = res['mcal_c']

        d['mcal_R'][i] = res['mcal_R']
        d['mcal_Rpsf'][i] = res['mcal_Rpsf']
        d['mcal_gpsf'][i] = res['mcal_gpsf']

    def _compare_psf_obs_fit(self, obs, **keys):
        gm = obs.get_gmix()
        d=gm._get_gmix_data()
        print("psf gmix F values:",d['p'])
        print("psf gmix T values:",d['irr']+d['icc'])
        print("psf gmix T rel:",(d['irr']+d['icc'])/gm.get_T())
        im=gm.make_image(obs.image.shape, jacobian=obs.jacobian)
        im *= obs.image.sum()/im.sum()

        self._compare_images(obs.image, im, **keys)

        key=raw_input("hit a key: ")
        if key=='q':
            stop

    def _compare_gal_obs_fit(self, obs, gm, **keys):
        d=gm._get_gmix_data()
        im=gm.make_image(obs.image.shape, jacobian=obs.jacobian)

        self._compare_images(obs.image, im, **keys)

        key=raw_input("hit a key: ")
        if key=='q':
            stop

 
    def _compare_images(self, im1, im2, **keys):
        import images
        keys['width']=1000
        keys['height']=1000
        images.compare_images(im1, im2, **keys)


class MaxMetacalSimnFitter(MaxMetacalFitter):

    def _do_metacal(self, boot):
        from ngmix import Bootstrapper

        super(MaxMetacalSimnFitter,self)._do_metacal(boot)
        print("    Calculating Rnoise")

        mb_obs_list = boot.mb_obs_list

        fitter = boot.get_max_fitter()
        res = fitter.get_result()

        gmix_list = [fitter.get_gmix()]

        # for noise added *before* metacal steps
        mobs_before = ngmix.simobs.simulate_obs(
            gmix_list,
            mb_obs_list,
            add_noise=True,
            convolve_psf=True
        )
        # for noise added *after* metacal steps
        mobs_after = ngmix.simobs.simulate_obs(
            gmix_list,
            mb_obs_list,
            add_noise=False,
            convolve_psf=True
        )

        boot_model_before=ngmix.Bootstrapper(mobs_before,
                                             use_logpars=self['use_logpars'],
                                             verbose=False)
        boot_model_after=ngmix.Bootstrapper(mobs_after,
                                            use_logpars=self['use_logpars'],
                                            verbose=False)


        mcal_obs_after = boot_model_after.get_metacal_obsdict(
            mobs_after[0][0],
            self['metacal_pars']
        )

        # now add noise after creating the metacal observations
        # using the same noise image!

        noise = mobs_before[0][0].noise_image
        for key in mcal_obs_after:
            obs=mcal_obs_after[key]
            obs.image = obs.image + noise

        super(MaxMetacalSimnFitter,self)._do_metacal(
            boot_model_before
        )
        super(MaxMetacalSimnFitter,self)._do_metacal(
            boot_model_after,
            metacal_obs=mcal_obs_after
        )

        res_before = boot_model_before.get_metacal_max_result()
        res_after = boot_model_after.get_metacal_max_result()

        Rnoise = res_before['mcal_R'] - res_after['mcal_R']
        Rpsf_noise = res_before['mcal_Rpsf'] - res_after['mcal_Rpsf']

        res['mcal_Rnoise'] = Rnoise
        res['mcal_Rpsf_noise'] = Rpsf_noise


    def _get_dtype(self):
        """
        get the dtype for the output struct
        """
        dt=super(MaxMetacalSimnFitter,self)._get_dtype()

        dt += [
            ('mcal_Rnoise','f8',(2,2)),
            ('mcal_Rpsf_noise','f8',2),
        ]
        return dt

    def _copy_to_output(self, res, i):
        """
        copy parameters specific to this class
        """
        super(MaxMetacalSimnFitter,self)._copy_to_output(res, i)

        d=self.data
        d['mcal_Rnoise'][i] = res['mcal_Rnoise']
        d['mcal_Rpsf_noise'][i] = res['mcal_Rpsf_noise']

 
class MaxMetacalFitterDegrade(MaxMetacalFitter):
    """
    degrade the images
    """

    def _setup(self, *args, **kw):
        super(MaxMetacalFitterDegrade,self)._setup(*args, **kw)

        if len(self['s2n_target']) > 1:
            raise ValueError("only one target s2n allowed for now")

        s2n_target = float(self['s2n_target'][0])
        self['noise_boost'] = self.sim['s2n_for_noise']/s2n_target
        self['extra_noise'] = self.sim['noise']*self['noise_boost']

        print("    boosting noise by",self['noise_boost'])

class MaxMetacalFitterDegradeGS(MaxMetacalFitterDegrade):
    """
    Add extra noise to the image
    """

    def _setup(self, *args, **kw):
        super(MaxMetacalFitterDegrade,self)._setup(*args, **kw)

        noise = self.sim['noise']

        extra_noise = sqrt(self['target_noise']**2 - noise**2)
        self['extra_noise'] = extra_noise

    def _do_metacal(self, boot):

        ppars=self['psf_pars']
        mconf=self['max_pars']
        Tguess=self.sim.get('psf_T',4.0)
        psf_fit_pars = ppars.get('fit_pars',None)

        extra_noise=self['extra_noise']
        nrand=self['nrand']
        print("    adding extra noise:",extra_noise, "nrand:",nrand)

        boot.fit_metacal_max_addnoise(
            extra_noise,
            nrand,
            ppars['model'],
            self['fit_model'],
            mconf['pars'],
            Tguess,
            psf_fit_pars=psf_fit_pars,
            prior=self.prior,
            ntry=mconf['ntry'],
            metacal_pars=self['metacal_pars'],
        )


class MaxMetanoiseFitter(MaxMetacalFitter):
    """
    degrade the images
    """

    def _do_fits(self, obs):
        """
        the basic fitter for this class
        """

        rdict=super(MaxMetanoiseFitter,self)._do_fits(obs)
        boot=rdict['boot']

        print("doing metanoise")

        ppars=self['psf_pars']
        mconf=self['max_pars']
        Tguess=self.sim.get('psf_T',4.0)
        psf_fit_pars = ppars.get('fit_pars',None)

        try:
            boot.fit_metanoise_max(ppars['model'],
                                   self['fit_model'],
                                   mconf['pars'],
                                   Tguess,
                                   self['metanoise_pars']['nrand'],
                                   psf_fit_pars=psf_fit_pars,
                                   prior=self.prior,
                                   ntry=mconf['ntry'],
                                   metacal_pars=self['metacal_pars'])
        except BootPSFFailure:
            raise TryAgainError("failed to fit metacal psfs")
        except BootGalFailure:
            raise TryAgainError("failed to fit metacal metanoise galaxy")

        fitter=boot.get_max_fitter() 
        res=fitter.get_result()

        mnres = boot.get_metanoise_max_result()
        res.update(mnres)

        return {'fitter':fitter,
                'boot':boot,
                'res':res}

    def _print_res(self,res):
        """
        print some stats
        """

        super(MaxMetanoiseFitter,self)._print_res(res)
        print("    mnoise s2n_r:",res['mnoise_s2n_r'],'mnoise s2n_simple:',res['mnoise_s2n_simple'])
        print_pars(res['mnoise_pars'],       front='    mnoise pars: ')
        print_pars(res['mnoise_R'].ravel(),  front='    mnoise R:    ')
        print_pars(res['mnoise_Rpsf'],       front='    mnoise Rpsf: ')

        print("    mnoise c:",res['mnoise_c'][0], res['mnoise_c'][1])


    def _get_dtype(self):
        """
        get the dtype for the output struct
        """
        dt=super(MaxMetanoiseFitter,self)._get_dtype()

        npars=self['npars']
        dt += [
            ('mnoise_pars','f8',npars),
            ('mnoise_pars_cov','f8',(npars,npars)),
            ('mnoise_g','f8',2),
            ('mnoise_g_cov','f8', (2,2) ),
            ('mnoise_g_noshear','f8',2),
            ('mnoise_c','f8',2),
            ('mnoise_s2n_r','f8'),
            ('mnoise_s2n_simple','f8'),
            ('mnoise_R','f8',(2,2)),
            ('mnoise_Rpsf','f8',2),
            ('mnoise_gpsf','f8',2),
        ]
        return dt


    def _copy_to_output(self, res, i):
        """
        copy parameters specific to this class
        """
        super(MaxMetanoiseFitter,self)._copy_to_output(res, i)

        d=self.data

        d['mnoise_pars'][i] = res['mnoise_pars']
        d['mnoise_pars_cov'][i] = res['mnoise_pars_cov']
        d['mnoise_g'][i] = res['mnoise_g']
        d['mnoise_g_cov'][i] = res['mnoise_g_cov']
        d['mnoise_s2n_r'][i] = res['mnoise_s2n_r']
        d['mnoise_s2n_simple'][i] = res['mnoise_s2n_simple']

        d['mnoise_g_noshear'][i] = res['mnoise_pars_noshear'][2:2+2]
        d['mnoise_c'][i] = res['mnoise_c']

        d['mnoise_R'][i] = res['mnoise_R']
        d['mnoise_Rpsf'][i] = res['mnoise_Rpsf']
        d['mnoise_gpsf'][i] = res['mnoise_gpsf']

def _make_new_obs(obs, im):
    """
    only diff is im
    """

    psf_obs=Observation(obs.psf.image,
                        jacobian=obs.psf.jacobian,
                        gmix=obs.psf.gmix.copy())

    newobs=Observation(im,
                       jacobian=obs.jacobian,
                       weight=obs.weight,
                       psf=psf_obs)
    return newobs




class EMMetacalFitter(SimpleFitterBase):
    """
    this is not finished
    """


    def _dofit(self, imdict):

        obs = imdict['obs']
        resdict = self._do_one_fit(obs)
        psf_resdict = self._do_one_fit(obs.psf)

        res=resdict['res']
        res['epsf'] = psf_resdict['res']['e']

        boot=ngmix.Bootstrapper(obs)

        mcal_pars = self['metacal_pars']
        mcal_obs_dict = boot.get_metacal_obsdict(obs, mcal_pars)

        epsf = zeros( (4,2) )
        ipsf=0

        pars_dict = {}
        for key in mcal_obs_dict:
            tobs = mcal_obs_dict[key]

            psfd = self._do_one_fit(tobs.psf, type='psf')
            rd = self._do_one_fit(tobs)

            if 'psf' not in key and 'noshear' not in key:
                epsf[ipsf,:] = psfd['res']['e']
                ipsf+=1

            pars_dict[key] = rd['res']['pars']


        Rres=self._calc_R(pars_dict)

        res.update(Rres)
        res['mcal_epsf'] = epsf.mean(axis=0)

        return resdict['fitter']

    def _calc_R(self, pars):

        step = self['metacal_pars']['step']

        pars_mean = (pars['1p']+
                     pars['1m']+
                     pars['2p']+
                     pars['2m'])/4.0


        R=zeros( (2,2) ) 
        Rpsf=zeros(2)

        fac = 1.0/(2.0*step)

        R[0,0] = (pars['1p'][2]-pars['1m'][2])*fac
        R[0,1] = (pars['1p'][3]-pars['1m'][3])*fac
        R[1,0] = (pars['2p'][2]-pars['2m'][2])*fac
        R[1,1] = (pars['2p'][3]-pars['2m'][3])*fac

        Rpsf[0] = (pars['1p_psf'][2]-pars['1m_psf'][2])*fac
        Rpsf[1] = (pars['2p_psf'][3]-pars['2m_psf'][3])*fac

        pars_noshear = pars['noshear']

        c = pars_mean[2:2+2] - pars_noshear[2:2+2]

        res = {'mcal_pars':pars_mean,
               'mcal_e':pars_mean[2:2+2],
               'mcal_R':R,
               'mcal_Rpsf':Rpsf,
               'mcal_e_noshear':pars_noshear[2:2+2].copy(),
               'mcal_c':c}
        return res


    def _do_one_fit(self, obs, type='gal'):
        """
        this is EM specific
        """
        from ngmix.bootstrap import EMRunner

        emconf = self['em_pars']
        ngauss = emconf['ngauss']

        ntry = emconf['ntry']
 
        Tguess=self._get_Tguess(ngauss, type)
        empars = emconf['pars']
        runner=EMRunner(obs, Tguess, ngauss, empars)

        runner.go(ntry=ntry)

        fitter=runner.get_fitter()
        res=fitter.get_result()

        if res['flags'] != 0:
            raise TryAgainError("em gal failed")

        res=fitter.get_result()
        res['model'] = self['fit_model']

        self._convert2pars(fitter)

        return {'boot':None,
                'fitter':fitter,
                'res':res}

    def _convert2pars(self, fitter):
        """
        convert em gauss to pars
        """
        gm=fitter.get_gmix()

        cen1,cen2=gm.get_cen()
        e1,e2,T = gm.get_e1e2T()

        pars=array([cen1,cen2, e1, e2, T, 1.0])

        res=fitter.get_result()
        res['pars']=pars
        res['e']=array([e1,e2])

    def _get_Tguess(self, ngauss, type):
        """
        Get the starting guess.
        """

        if type=='gal':
            Tguess=6.0
        else:
            Tguess=4.0
        return Tguess*(1.0 + 0.05*srandu())

    def _set_prior(self):
        """
        Set all the priors
        """
        pass
 
    def _print_res(self,res):
        """
        print some stats
        """

        mess="    ntry: %d  numiter: %d"
        mess = mess % (res['ntry'],res['numiter'])
        print(mess)

        print_pars(res['pars'],            front='    pars: ')

        print_pars(res['mcal_pars'],       front='    mcal pars: ')
        print_pars(res['mcal_R'].ravel(),  front='    mcal R:    ')
        print_pars(res['mcal_Rpsf'],       front='    mcal Rpsf: ')

        print("    mcal c:",res['mcal_c'][0], res['mcal_c'][1])



    def _copy_to_output(self, res, i):
        """
        this is EM specific due to numiter, fdiff
        """

        super(SimpleFitterBase,self)._copy_to_output(res, i)

        d=self.data

        for f in ['pars','e','epsf','numiter','fdiff']:
            d[f][i] = res[f]

        for f in ['pars','e','epsf','e_noshear',
                  'R','Rpsf','c']:

            f = 'mcal_%s' % f
            d[f][i] = res[f]

        

    def _get_dtype(self):
        """
        this is EM specific due to numiter, fdiff
        """
        dt=super(SimpleFitterBase,self)._get_dtype()

        npars=self['npars']
        dt += [
            ('pars','f8',npars),
            ('e','f8',2),
            ('epsf','f8',2),
            ('numiter','i4'),
            ('fdiff','f8'),
            ('mcal_pars','f8',npars),
            ('mcal_e','f8',2),
            ('mcal_epsf','f8',2),
            ('mcal_e_noshear','f8',2),
            ('mcal_R','f8',(2,2)),
            ('mcal_Rpsf','f8',2),
            ('mcal_c','f8',2)
        ]
        return dt



class NCalFitter(MaxFitter):
    def _dofit(self, imdict):
        """
        Fit according to the requested method
        """

        self.imdict=imdict


        obs=imdict['obs']
        mdict = self._do_fits(obs)


        res=mdict['res']

        nres=self._do_ncal(mdict['boot'])

        res.update(nres)

        fitter=mdict['fitter']
        self.fitter=fitter

        return fitter

    def _do_fits(self, obs, fit_psfs=True):
        """
        the basic fitter for this class
        """
        from ngmix import Bootstrapper

        boot=Bootstrapper(obs,
                          use_logpars=self['use_logpars'],
                          verbose=False)

        if fit_psfs:
            Tguess=self.sim.get('psf_T',4.0)
            ppars=self['psf_pars']

            psf_fit_pars = ppars.get('fit_pars',None)
        
            try:
                boot.fit_psfs(ppars['model'],
                              Tguess,
                              ntry=ppars['ntry'],
                              fit_pars=psf_fit_pars)
            except BootPSFFailure:
                raise TryAgainError("failed to fit psf")

        mconf=self['max_pars']

        try:
            boot.fit_max(self['fit_model'],
                         mconf['pars'],
                         prior=self.prior,
                         ntry=mconf['ntry'])


            if mconf['pars']['method']=='lm':
                boot.try_replace_cov(mconf['cov_pars'])


        except BootPSFFailure:
            raise TryAgainError("failed to fit ncal psfs")
        except BootGalFailure:
            raise TryAgainError("failed to fit ncal galaxy")

        fitter=boot.get_max_fitter() 
        res=fitter.get_result()

        return {'fitter':fitter,
                'boot':boot,
                'res':res}

    def _get_simulated_obs(self, boot):

        oobs = boot.mb_obs_list[0][0]

        mfitter=boot.get_max_fitter()
        model = mfitter.model_name

        mres=mfitter.get_result()
        mpars = mres['pars'].copy()

        if self['use_logpars']:
            mpars[4] = exp(mpars[4])
            mpars[5] = exp(mpars[5])

        g1=mpars[2]
        g2=mpars[3]
        T0=mpars[4]

        T0round = ngmix.moments.get_Tround(T0, g1, g2)

        sh0 = Shape(g1,g2)

        step = self['ncal_pars']['step']
        typelist = [
            ('1p', Shape( step,  0.0)),
            ('1m', Shape(-step,  0.0)),
            ('2p', Shape( 0.0,   step)),
            ('2m', Shape( 0.0,  -step))
        ]

        noise_image = ngmix.simobs.get_noise_image(oobs.weight)

        ndict={}
        for i,t in enumerate(typelist):
            type=t[0]
            shear=t[1]

            sh = sh0.get_sheared(shear) 
            T  = ngmix.moments.get_T(T0round, sh.g1, sh.g2)

            pars = mpars.copy()
            pars[2] = sh.g1
            pars[3] = sh.g2
            pars[4] = T

            gm = GMixModel(pars, model)
            
            sobs = ngmix.simobs.simulate_obs(gm, oobs, add_noise=False)

            sobs.image = sobs.image + noise_image

            ndict[type] = sobs

        if False:
            import images
            images.compare_images(oobs.image, sobs.image,
                                  label1='image',
                                  label2='sim')

            key=raw_input('hit a key: ')
            if key=='q':
                stop

        return ndict

    def _do_ncal(self, boot):

        ndict = self._get_simulated_obs(boot)

        pdict={}

        for key in ndict:
            tres = self._do_fits(ndict[key], fit_psfs=False)

            tpars = tres['res']['pars']

            pdict[key] = tpars

        nres={}

        R=zeros( (2,2) ) 
        Rpsf=zeros(2)

        step = self['ncal_pars']['step']
        fac = 1.0/(2.0*step)

        R[0,0] = (pdict['1p'][2]-pdict['1m'][2])*fac
        R[0,1] = (pdict['1p'][3]-pdict['1m'][3])*fac
        R[1,0] = (pdict['2p'][2]-pdict['2m'][2])*fac
        R[1,1] = (pdict['2p'][3]-pdict['2m'][3])*fac

        nres['ncal_R'] = R
        nres['ncal_Rpsf'] = Rpsf

        return nres

    def _print_res(self,res):
        """
        print some stats
        """

        super(NCalFitter,self)._print_res(res)
        print_pars(res['ncal_R'].ravel(),  front='    ncal R:    ')

    def _get_dtype(self):
        """
        get the dtype for the output struct
        """
        dt=super(NCalFitter,self)._get_dtype()

        npars=self['npars']
        dt += [
            ('ncal_R','f8',(2,2)),
            ('ncal_Rpsf','f8',2),
        ]
        return dt


    def _copy_to_output(self, res, i):
        """
        copy parameters specific to this class
        """
        super(NCalFitter,self)._copy_to_output(res, i)

        d=self.data

        d['ncal_R'][i] = res['ncal_R']
        d['ncal_Rpsf'][i] = res['ncal_Rpsf']



def make_sheared_pars(pars, shear_g1, shear_g2):
    from ngmix import Shape
    shpars=pars.copy()

    sh=Shape(shpars[2], shpars[3])
    sh.shear(shear_g1, shear_g2)

    shpars[2]=sh.g1
    shpars[3]=sh.g2

    return shpars

def make_cheat_metacal_obs(psf_gmix, pars, model, noise_image, obs, nsub):
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


