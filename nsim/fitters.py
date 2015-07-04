"""
fit simulated images
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
from ngmix import srandu
from ngmix.fitting import print_pars
from ngmix.gexceptions import GMixRangeError
from ngmix.observation import Observation
from ngmix.gexceptions import GMixMaxIterEM
from ngmix.shape import Shape

from ngmix.gexceptions import BootPSFFailure, BootGalFailure

from .sim import NGMixSim
from .util import write_fits

# minutes
DEFAULT_CHECKPOINTS=[30,60,90,110]

class TryAgainError(Exception):
    def __init__(self, message):

        # Call the base class constructor with the parameters it needs
        Exception.__init__(self, message)


class FitterBase(dict):
    def __init__(self, sim, run_conf, npairs, **keys):

        self.sim=sim
        self._setup(run_conf, **keys)

        self['npairs']=npairs

        self._set_prior()

        self._setup_checkpoints(**keys)
        if self.data is None:
            self._make_struct()

        pprint.pprint(self)

    def go(self):
        """
        process the requested number of pairs
        """

        self._start_timer()

        npairs=self['npairs']
        i=0
        for ipair in xrange(npairs):
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
                        print(str(err))

                self._copy_to_output(reslist[0], i)
                i += 1
                self._copy_to_output(reslist[1], i)
                i += 1

            self._set_elapsed_time()
            self._try_checkpoint()

        self._set_elapsed_time()

        print('time minutes:',self.tm_minutes)
        print('time per pair sec:',self.tm/npairs)
        print('time per image sec:',self.tm/(2*npairs))

    def process_pair(self):
        """
        Create a simulated image pair and perform the fit
        """
        imdicts = self.sim.get_image_pair()

        print(imdicts['im1']['obs'].image.shape)

        reslist=[]
        for key in ['im1','im2']:

            imd = imdicts[key]

            fitter, psf_flux_res=self._dofit(imd)

            res=fitter.get_result()
            if res['flags'] != 0:
                raise TryAgainError("failed at %s "
                                    "flags %s" % (key,res['flags']))

            res['pars_true'] = imd['pars']
            res['model_true'] = imd['model']
            res['s2n_true'] = imd['s2n']

            res['psf_flux'] = psf_flux_res['psf_flux']
            res['psf_flux_err'] = psf_flux_res['psf_flux_err']
            res['psf_flux_s2n'] = res['psf_flux']/res['psf_flux_err']

            self._print_res(res)

            reslist.append(res)

            if self['make_plots']:
                self._make_plots(fitter,key)

        return reslist


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

    def _setup(self, run_conf, **keys):
        """
        Check and set the configurations
        """


        if self.sim['name'] != run_conf['sim']:
            err="sim name in run config '%s' doesn't match sim name '%s'"
            raise ValueError(err % (run_conf['sim'],self.sim['name']))

        self.update(run_conf)

        self['use_logpars']=self.get('use_logpars',False)
        self['npars']=ngmix.gmix.get_model_npars(self['fit_model'])
        self['npars_true'] = ngmix.gmix.get_model_npars(self.sim['obj_model'])

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
           ]

        return dt

    def _make_struct(self):
        """
        Make the output array
        """

        dt=self._get_dtype()
        self.data=numpy.zeros(self['npairs']*2, dtype=dt)

class SimpleFitterBase(FitterBase):

    def _copy_to_output(self, res, i):

        super(SimpleFitterBase,self)._copy_to_output(res, i)

        d=self.data

        d['psf_flux'][i] = res['psf_flux']
        d['psf_flux_err'][i] = res['psf_flux_err']
        d['psf_flux_s2n'][i] = res['psf_flux_s2n']

        d['pars'][i,:] = res['pars']
        d['pcov'][i,:,:] = res['pars_cov']

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
        else:
            T_prior_pars = Tp['pars']
            fit_T_prior=ngmix.priors.TwoSidedErf(*T_prior_pars)

        cp=ppars['counts']
        if cp['type']=="truth":
            print("using true counts pdf for prior")

            if self['use_logpars']:

                print("    converting to log")
                counts       = self.sim['obj_counts_mean']
                counts_sigma = self.sim['obj_counts_sigma_frac']*counts
                logc_mean, logc_sigma=ngmix.priors.lognorm_convert(counts,
                                                                   counts_sigma)
                fit_counts_prior = ngmix.priors.Normal(logc_mean, logc_sigma)

            else:
                fit_counts_prior = self.sim.counts_pdf

        else:
            counts_prior_pars = cp['pars']
            fit_counts_prior=ngmix.priors.TwoSidedErf(*counts_prior_pars)

        cp=ppars['cen']
        if cp['type']=="truth":
            fit_cen_prior=self.sim.cen_pdf
        else:
            fit_cen_sigma_arcsec=cp['sigma']
            fit_cen_prior=ngmix.priors.CenPrior(0.0,
                                                0.0,
                                                fit_cen_sigma_arcsec,
                                                fit_cen_sigma_arcsec)


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
            ('psf_flux','f8'),
            ('psf_flux_err','f8'),
            ('psf_flux_s2n','f8'),
            ('pars','f8',npars),
            ('pcov','f8',(npars,npars)),
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

        boot=ngmix.Bootstrapper(imdict['obs'], use_logpars=self['use_logpars'])

        Tguess=self.sim['psf_T']
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
        psf_flux_res = boot.get_psf_flux_result()

        return fitter, psf_flux_res

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

        mpars=self['metacal_pars']

        obs=imdict['obs']
        mdict = self._do_one_fit(obs)
        res=mdict['res']
        psf_flux_res=mdict['psf_flux_res']
        fitter=mdict['fitter']

        self.fitter=fitter

        sdict = self._get_sensitivity(obs)
        res.update(sdict)

        g_sens_model = self._get_sensitivity_model(obs, fitter)
        res['g_sens_model'] = g_sens_model

        return fitter, psf_flux_res

    def _do_one_fit(self, obs, guess=None, ntry=None):
        """
        the basic fitter for this class
        """
        boot=ngmix.Bootstrapper(obs, use_logpars=self['use_logpars'])

        Tguess=self.sim['psf_T']
        ppars=self['psf_pars']
        try:
            boot.fit_psfs(ppars['model'], Tguess, ntry=ppars['ntry'])
        except BootPSFFailure:
            raise TryAgainError("failed to fit psf")

        mconf=self['max_pars']
        if ntry is None:
            ntry=mconf['ntry']
        try:
            boot.fit_max(self['fit_model'],
                         mconf['pars'],
                         prior=self.prior,
                         guess=guess,
                         ntry=mconf['ntry'])

            if mconf['pars']['method']=='lm':
                boot.try_replace_cov(mconf['cov_pars'])
        except BootGalFailure:
            raise TryAgainError("failed to fit galaxy")

        fitter=boot.get_max_fitter() 
        res=fitter.get_result()
        psf_flux_res = boot.get_psf_flux_result()

        return {'fitter':fitter,
                'boot':boot,
                'res':res,
                'psf_flux_res':psf_flux_res}

    def _get_sensitivity(self, obs):
        """
        fit sheared observations
        """
        mpars=self['metacal_pars']

        R_obs_1m, R_obs_1p, R_obs_noshear = \
                self._get_metacal_obslist(obs, get_noshear=True)

        mdict_noshear = self._do_one_fit(R_obs_noshear)
        res_noshear = mdict_noshear['res']
        pars_noshear=res_noshear['pars'].copy()


        mdict_1m = self._do_one_fit(R_obs_1m, guess=pars_noshear, ntry=1)
        mdict_1p = self._do_one_fit(R_obs_1p, guess=pars_noshear, ntry=1)
        res_1m = mdict_1m['res']
        res_1p = mdict_1p['res']

        pars_1m=res_1m['pars']
        pars_1p=res_1p['pars']
        pars_mean = 0.5*(pars_1m + pars_1p)
        
        g_1m=res_1m['pars'][2:2+2]
        g_1p=res_1p['pars'][2:2+2]

        step=mpars['step']
        g_sens1 = (g_1p[0] - g_1m[0])/(2*step)

        g_mean = pars_mean[2:2+2]
        c = g_mean-pars_noshear[2:2+2]

        print("    c: %g %g" % tuple(c))

        g_sens = array( [g_sens1]*2 )


        return {'pars_noshear':pars_noshear,
                'pars_mean':pars_mean,
                'g_sens':g_sens,
                'c':c}

    def _get_sensitivity_model(self, obs, fitter):
        """
        just simulate the model at different shears

        would only correct noise bias I think
        """
        mpars=self['metacal_pars']

        im=obs.image
        imshape=im.shape
        jacobian=obs.jacobian

        res=fitter.get_result()
        model=res['model']

        pars=fitter.get_band_pars(res['pars'], band=0)
        #gm0=ngmix.GMixModel(pars,model)
        #gm=gm0.convolve(obs.psf.gmix)
        #model_im=gm.make_image(imshape, jacobian=jacobian, nsub=1)
        #nim=obs.image-model_im

        skysig=self.sim['skysig']
        nim=numpy.random.normal(scale=skysig, size=imshape)

        parsm=pars.copy()
        parsp=pars.copy()

        shm=ngmix.Shape(parsm[2], parsm[3])
        shp=ngmix.Shape(parsp[2], parsp[3])
        
        step=mpars['step']
        shm.shear(-step, 0.0)
        shp.shear( step, 0.0)

        parsm[2],parsm[3]=shm.g1, shm.g2
        parsp[2],parsp[3]=shp.g1, shp.g2

        gmm0=ngmix.GMixModel(parsm,model)
        gmp0=ngmix.GMixModel(parsp,model)

        gmm=gmm0.convolve(obs.psf.gmix)
        gmp=gmp0.convolve(obs.psf.gmix)

        imm=gmm.make_image(imshape, jacobian=jacobian, nsub=1)
        imp=gmp.make_image(imshape, jacobian=jacobian, nsub=1)

        imm += nim
        imp += nim
        
        R_obs_1m = _make_new_obs(obs, imm)
        R_obs_1p = _make_new_obs(obs, imp)

        mdict_1m = self._do_one_fit(R_obs_1m)
        mdict_1p = self._do_one_fit(R_obs_1p)
        res_1m = mdict_1m['res']
        res_1p = mdict_1p['res']

        pars_1m=res_1m['pars']
        pars_1p=res_1p['pars']
        g_1m=pars_1m[2:2+2]
        g_1p=pars_1p[2:2+2]

        g_sens1 = (g_1p[0] - g_1m[0])/(2*step)

        g_sens = array( [g_sens1]*2 )

        return g_sens


    def _get_metacal_obslist(self, obs, get_noshear=False):
        """
        get Observations for the sheared images
        """
        from ngmix.metacal import Metacal

        mpars=self['metacal_pars']

        mc=Metacal(obs,
                   whiten=mpars['whiten'],
                   same_seed=mpars['same_seed'])

        sval=mpars['step']
        sh1m=Shape(-sval,  0.00 )
        sh1p=Shape( sval,  0.00 )
        #sh2m=Shape( 0.00, -sval )
        #sh2p=Shape( 0.00,  sval )


        R_obs1p = mc.get_obs_galshear(sh1p)
        if get_noshear:
            R_obs1m, R_obs_noshear = mc.get_obs_galshear(sh1m, get_unsheared=True)
            return R_obs1m, R_obs1p, R_obs_noshear 
        else:
            R_obs1m = mc.get_obs_galshear(sh1m)
            return R_obs1m, R_obs1p

    def _print_res(self,res):
        """
        print some stats
        """

        super(MaxMetacalFitter,self)._print_res(res)
        print_pars(res['g_sens'],      front='        sens: ')

    def _setup(self, run_conf, **keys):
        """
        setup parameters for this specific class
        """
        super(MaxMetacalFitter,self)._setup(run_conf, **keys)
        mpars=self['metacal_pars']
        mpars['whiten']=mpars.get('whiten',False)
        mpars['same_seed']=mpars.get('same_seed',False)


    def _get_dtype(self):
        """
        get the dtype for the output struct
        """
        dt=super(MaxMetacalFitter,self)._get_dtype()
        dt += [
            ('g_noshear','f8',2),
            ('g_mean','f8',2),
            ('g_sens','f8',2),
            ('g_sens_model','f8',2),
        ]
        return dt


    def _copy_to_output(self, res, i):
        """
        copy parameters specific to this class
        """
        super(MaxMetacalFitter,self)._copy_to_output(res, i)

        d=self.data
        # reconv with no shear
        d['g_noshear'][i] = res['pars_noshear'][2:2+2]

        # mean from the two sheared images
        d['g_mean'][i] = res['pars_mean'][2:2+2]

        # sensitivity from two sheared images
        d['g_sens'][i] = res['g_sens']

        # sensitivity from simulating the model
        d['g_sens_model'][i] = res['g_sens_model']


class AdmomMetacalFitter(MaxMetacalFitter):
    def _do_one_fit(self, obs, **kw):
        """
        note guess= and ntry= are ignored
        """
        import admom
        admompars=self['admom_pars']
        ntry=admompars['ntry']

        gm=self.imdict['gm']

        row0,col0=obs.jacobian.get_cen()
        row,col=gm.get_cen()
        row += row0
        col += col0

        T=gm.get_T()

        for i in xrange(ntry):
            rowguess=row + 0.1*srandu()
            colguess=col + 0.1*srandu()
            Tguess=T*(1.0 + 0.1*srandu())

            res=admom.admom(obs.image, rowguess, colguess,
                            sigsky=self.sim['skysig'],
                            guess=Tguess/2.,
                            **admompars)

            res['flags']=res['whyflag']
            if res['flags']==0:
                break

        if res['flags'] != 0:
            raise TryAgainError("admom error '%s'" % res['whystr'])

        res['pars'] = array([res['wrow']-row0,
                             res['wcol']-col0,
                             res['e1'],
                             res['e2'],
                             res['Irr']+res['Icc'],
                             1.0])

        res['pars_err']=res['pars']*0 - 9999

        res['pars_cov'] = diag( res['pars']*0 - 9999 )

        res['g'] = res['pars'][2:2+2].copy()
        res['g_cov'] = res['pars_cov'][2:2+2, 2:2+2].copy()

        res['s2n_w'] = res['s2n']

        psf_flux_res={'psf_flux':-9999.0,
                      'psf_flux_err':9999.0}

        reswrap=MomWrapper(res)
        return {'fitter':reswrap,
                'res':res,
                'psf_flux_res':psf_flux_res}

    def _get_sensitivity_model(self, obs, fitter):
        return zeros(2)-9999.0

    #def _copy_to_output(self, res, i):
    #    super(SimpleFitterBase,self)._copy_to_output(res, i)


class MomMetacalFitter(MaxMetacalFitter):
    def _do_one_fit(self, obs, ntry=None, **kw):
        """
        note guess= is ignored
        """

        mompars=self['mom_pars']

        wtype=mompars['weight_type']
        if wtype =="truth":
            gm_weight = self.imdict['gm'].copy()
            # weight should always have center 0,0, coinciding
            # with jacobian center of coordinates.
            gm_weight.set_cen(0.0, 0.0)
        else:
            raise RuntimeError("bad weight type: '%s'" % wtype)

        res=gm_weight.get_weighted_mom_sums(obs, **mompars)

        if res['flags'] != 0:
            raise TryAgainError("error getting weights: '%s'" % res['flagstr'])

        skyvar=self.sim['skysig']**2

        #res['pars'][2] = res['pars'][2]/res['pars'][5]
        #res['pars'][3] = res['pars'][3]/res['pars'][5]
        res['pars'][2] = res['pars'][2]/res['pars'][4]
        res['pars'][3] = res['pars'][3]/res['pars'][4]

        res['pars_var'] *= skyvar
        res['pars_cov'] *= skyvar
        res['pars_err'] = sqrt(res['pars_var'])
        res['g'] = res['pars'][2:2+2].copy()
        res['g_cov'] = res['pars_cov'][2:2+2, 2:2+2].copy()

        res['s2n_w'] = res['pars'][5]/res['pars_err'][5]

        psf_flux_res={'psf_flux':-9999.0,
                      'psf_flux_err':9999.0}

        reswrap=MomWrapper(res)
        return {'fitter':reswrap,
                'res':res,
                'psf_flux_res':psf_flux_res}

    def _get_sensitivity_model(self, obs, fitter):
        return zeros(2)-9999.0

    #def _copy_to_output(self, res, i):
    #    super(SimpleFitterBase,self)._copy_to_output(res, i)


class MomWrapper(object):
    def __init__(self, res):
        self._result=res
    def get_result(self):
        return self._result

class MaxMetacalFitterModel(MaxMetacalFitter):
    def __init__(self, *args, **keys):
        super(MaxMetacalFitterModel,self).__init__(*args,**keys)
        mpars=self['metacal_pars']
        if mpars['mean_from'] != 'orig':
            raise ValueError("mean_from must be orig")

    def _get_sensitivity(self, obs):
        """
        just simulate the model at different shears

        would only correct noise bias I think
        """
        mpars=self['metacal_pars']

        im=obs.image
        imshape=im.shape
        jacobian=obs.jacobian

        res=self.fitter.get_result()
        model=res['model']

        pars=self.fitter.get_band_pars(res['pars'], band=0)
        #gm0=ngmix.GMixModel(pars,model)
        #gm=gm0.convolve(obs.psf.gmix)
        #model_im=gm.make_image(imshape, jacobian=jacobian, nsub=1)
        #nim=obs.image-model_im

        skysig=self.sim['skysig']
        nim=numpy.random.normal(scale=skysig, size=imshape)

        parsm=pars.copy()
        parsp=pars.copy()

        shm=ngmix.Shape(parsm[2], parsm[3])
        shp=ngmix.Shape(parsp[2], parsp[3])
        
        step=mpars['step']
        shm.shear(-step, 0.0)
        shp.shear( step, 0.0)

        parsm[2],parsm[3]=shm.g1, shm.g2
        parsp[2],parsp[3]=shp.g1, shp.g2

        gmm0=ngmix.GMixModel(parsm,model)
        gmp0=ngmix.GMixModel(parsp,model)

        gmm=gmm0.convolve(obs.psf.gmix)
        gmp=gmp0.convolve(obs.psf.gmix)

        imm=gmm.make_image(imshape, jacobian=jacobian, nsub=1)
        imp=gmp.make_image(imshape, jacobian=jacobian, nsub=1)


        imm += nim
        imp += nim
        
        R_obs_1m = _make_new_obs(obs, imm)
        R_obs_1p = _make_new_obs(obs, imp)

        mdict_1m = self._do_one_fit(R_obs_1m)
        mdict_1p = self._do_one_fit(R_obs_1p)
        res_1m = mdict_1m['res']
        res_1p = mdict_1p['res']

        pars_1m=res_1m['pars']
        pars_1p=res_1p['pars']
        g_1m=pars_1m[2:2+2]
        g_1p=pars_1p[2:2+2]

        g_sens1 = (g_1p[0] - g_1m[0])/(2*step)

        g_sens = array( [g_sens1]*2 )

        return None, g_sens

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





class ISampleMetacalFitter(MaxMetacalFitter):
    def _dofit(self, imdict):
        """
        Fit according to the requested method
        """
        mpars=self['metacal_pars']

        obs=imdict['obs']

        # we want to create a sampler and covariance matrix
        # then fit the +/- shears using the same samples
        print("    Doing main fit")
        res, sampler, psf_flux_res = self._do_main_fit(obs)

        assert mpars['mean_from'] in ['orig','avg'],"avg or orig for now"
        print("    getting sens")
        pars, g_sens = self._get_sensitivity_avg(obs)


        if mpars['mean_from'] != 'orig':
            res['pars'] = pars
            res['g'] = pars[2:2+2].copy()
        else:
            print("    using orig g")
        res['g_sens'] = g_sens

        # reset result
        sampler._result=res
        return sampler, psf_flux_res

    def _do_main_fit(self, obs, ntry=None):
        mconf=self['max_pars']
        if ntry is None:
            ntry=mconf['ntry']

        # this does the max like fit only
        mdict = super(ISampleMetacalFitter,self)._do_one_fit(obs, ntry=ntry)
        res=mdict['res']
        psf_flux_res = mdict['psf_flux_res']
        boot= mdict['boot']

        ipars=self['isample_pars']

        try:
            boot.isample(ipars, prior=self.prior)
        except BootGalFailure:
            raise TryAgainError("failed to isample galaxy")

        self.boot=boot

        sampler=boot.get_isampler()

        maxres=boot.get_max_fitter().get_result()
        res['model'] = maxres['model']
        res['s2n_w'] = maxres['s2n_w']


        return res, sampler, psf_flux_res

    def _do_one_fit(self, obs, ntry=None):
        """
        re-use sampler and samples

        use the max fitter with new observation set
        """

        mconf=self['max_pars']
        if ntry is None:
            ntry=mconf['ntry']

        # just want to fit the psfs and get a maxlike fitter object
        mdict= super(ISampleMetacalFitter,self)._do_one_fit(obs, ntry=ntry)

        tboot=mdict['boot']
        print("    re-using samples")

        max_fitter=tboot.get_max_fitter()

        sampler=self.boot.get_isampler()
        sampler.set_iweights(max_fitter.calc_lnprob)
        sampler.calc_result()

        res=sampler.get_result()

        return {'res':res}

class EMMetacalFitter(MaxMetacalFitter):

    def _do_one_fit(self, obs, guess=None, ntry=None):
        from ngmix.em import GMixEM, prep_image
        from ngmix.bootstrap import EMRunner

        boot=None
        psf_flux_res={'psf_flux':-9999.0,
                      'psf_flux_err':9999.0}

        #if guess is not None:
        #    raise RuntimeError("support geuss")
        #    mdict = self._do_one_fit_with_guess(obs, guess)
        #    fitter=mdict['fitter']
        #    if ok:
        #        return boot, fitter, psf_flux_res

        emconf = self['em_pars']
        ngauss = emconf['ngauss']

        if ntry is None:
            ntry = emconf['ntry']
 
        Tguess=self._get_Tguess(ngauss)
        empars = emconf['pars']
        runner=EMRunner(obs, Tguess, ngauss, empars)

        runner.go(ntry=ntry)

        fitter=runner.get_fitter()
        res=fitter.get_result()
        res['model'] = 'gauss'

        if res['flags'] != 0:
            raise TryAgainError("em gal failed")

        self._convert2pars(fitter)

        return {'boot':boot,
                'fitter':fitter,
                'res':res,
                'psf_flux_res':psf_flux_res}

    def _do_one_fit_with_guess(self, obs, guess):
        from ngmix.em import GMixEM, prep_image

        image=obs.image
        imsky,sky=prep_image(image)

        sky_obs = Observation(imsky, jacobian=obs.jacobian)

        fitter=GMixEM(sky_obs)

        empars=self['em_pars']['pars']
        tol=empars['tol']
        maxiter=empars['maxiter']
 
        fitter.go(guess, sky,
                  tol=empars['tol'],
                  maxiter=empars['maxiter'])

        res=fitter.get_result()
        res['ntry'] = 1
        if res['flags']==0:
            self._convert2pars(fitter)
            ok=True
        else:
            ok=False

        return fitter, ok

    def _convert2pars(self, fitter):
        """
        convert em gauss to pars
        """
        gm=fitter.get_gmix()
        pars=gm.get_full_pars()

        T=pars[3]+pars[5]
        e1 = (pars[5]-pars[3])/T
        e2 = 2*pars[4]/T
        #g1,g2=e1,e2
        try:
            g1,g2=ngmix.shape.e1e2_to_g1g2(e1,e2)
        except GMixRangeError:
            raise TryAgainError("bad e1,e2")

        tpars=array([pars[1],pars[2], g1, g2, T, 1.0])

        res=fitter.get_result()
        res['g']=array([g1,g2])
        res['g_cov']=diag([9999.0,9999.0])
        res['pars']=tpars
        res['pars_err']=pars*0 + 9999.0
        res['pars_cov']=diag(pars*0 + 9999)
        res['s2n_w']=-9999.0
        res['nfev']=res['numiter']

    def _get_Tguess(self, ngauss):
        """
        Get the starting guess.
        """
        from .sim import srandu

        pt=self.sim.psf_gmix_true
        ngauss_true=len(pt)
        if ngauss is None or ngauss==ngauss_true:
            # we can just use the "truth" as a guess
            Tguess=pt.get_T()*(1.0 + 0.1*srandu())
        else:
            # more generic guess
            Tguess = pt.get_T()

            if ngauss==2:
                p1 = 0.6*(1.0 + 0.05*srandu() )
                p2 = 0.4*(1.0 + 0.05*srandu() )
                T1 = T*0.4*(1.0 + 0.1*srandu() )
                T2 = T*0.6*(1.0 + 0.1*srandu() )

                Tguess = (p1*T1 +p2*T2)/(p1+p2)

            else: 
                raise ValueError("support ngauss > 2")
            
        return Tguess


    def _get_sensitivity_model(self, obs, fitter):
        return zeros(2)-9999.0
