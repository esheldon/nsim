"""
fit simulated images
"""
from __future__ import print_function
try:
    xrange
except:
    xrange=range
    raw_input=input

import logging
import os
import time
import pprint

import numpy
from numpy import array, zeros, ones, log, log10, exp, sqrt, diag
from numpy import where, isfinite

import ngmix
from .util import log_pars
from ngmix.gexceptions import GMixRangeError
from ngmix.observation import Observation
from ngmix.gexceptions import GMixMaxIterEM
from ngmix.shape import Shape
from ngmix.gmix import GMixModel

from ngmix.gexceptions import BootPSFFailure, BootGalFailure

from .sim import NGMixSim

from . import util
from .util import write_fits, TryAgainError, load_gmixnd
from . import files

from copy import deepcopy

import esutil as eu

from . import deblending

try:
    import galsim
except ImportError:
    pass

logger = logging.getLogger(__name__)

class FitterBase(dict):
    def __init__(self, sim, run_conf, ngal, **keys):

        self.sim=sim
        self.rng=numpy.random.RandomState(
            seed=sim.rng.randint(0,2**30),
        )

        self._setup(run_conf, **keys)

        self['ngal']=ngal
        self['use_round_T'] = self.get('use_round_T',False)

        self.prior = self._get_prior()
        if 'masking' in self:
            self.replace_prior = self._get_prior(self['masking']['priors'])

        logger.info( pprint.pformat(self) )

    def go(self):
        """
        process the requested number of simulated galaxies
        """

        dlist = []
        self._start_timer()
        self.tm_sim=0.0
        self.tm_fit=0.0

        # total we attempted to fit
        n_proc=0

        detconf=self.get('detection',None)

        for igal in xrange(self['ngal']):
            self.igal=igal
            logger.debug('%s/%s' % (igal+1,self['ngal']) )

            try:

                tm0=time.time()
                obslist = self.sim()
                self.tm_sim += time.time()-tm0

                if self['show']:
                    import images
                    images.multiview(obslist[0].image)
                    if raw_input('hit a key (q to quit): ')=='q':
                        stop

                n_proc += 1

                # only objects that successfully get fit
                # contribute to this
                if detconf is not None and detconf['detect']:
                    tm0=time.time()
                    list_of_obslist=self._detect_objects(obslist)
                    this_dlist=self.process_multiple_objects(list_of_obslist)
                    dlist += this_dlist
                    self.tm_fit += time.time()-tm0
                else:
                    tm0=time.time()
                    res=self.process_one(obslist)
                    self.tm_fit += time.time()-tm0

                    d = self._make_output(res)

                    dlist.append(d)

                self._set_elapsed_time()

            except TryAgainError as err:
                logger.debug(str(err))

        if len(dlist) == 0:
            raise RuntimeError("none succeeded")

        self.data = eu.numpy_util.combine_arrlist(dlist)
        self._set_elapsed_time()

        logger.info("nprocessed: %s" % n_proc)
        logger.info('time minutes: %s' % self.tm_minutes)
        logger.info('time per (total) %s' % (self.tm/self['ngal']))
        logger.info('time to simulate: %s' % (self.tm_sim/self['ngal']))
        logger.info('time to fit: %s' % (self.tm_fit/n_proc))

    def process_one(self, obslist):
        """
        perform the fit
        """

        fitter=self._dofit(obslist)

        if isinstance(fitter, dict):
            res=fitter
        else:
            res=fitter.get_result()

        if res['flags'] != 0:
            raise TryAgainError("failed with flags %s" % res['flags'])

        self._print_res(res)

        if self['make_plots']:
            self._make_plots(fitter,self.get('fit_model',''))

        return res

    def process_multiple_objects(self, list_of_obslist):
        """
        process multiple objects independently
        """
        dlist=[]
        for obslist in list_of_obslist:
            try:
                res=self.process_one(obslist)
                d = self._make_output(res)

                dlist.append(d)
            except TryAgainError as err:
                logger.info(str(err))

        return dlist

    def get_data(self):
        """
        Get a ref to the data array with the fit results
        """
        return self.data

    def _make_plots(self, fitter, key):
        """
        Write a plot file of the trials
        """
        return

    def _dofit(self, imdict):
        """
        Fit according to the requested method
        """

        raise RuntimeError("over-ride me")


    def _detect_objects(self, obslist):
        from . import runsep

        assert len(obslist)==1
        obs=obslist[0]

        detconf=self['detection']

        objs, seg = runsep.find_objects(
            obs,
            segmentation_map=True,
        )

        if self['show']:
            self._show_image_and_seg(obs.image, seg)

        if objs.size == 0:
            logger.debug('no objects found')
            raise TryAgainError('no objects found')

        logger.debug('    found %d objects' % objs.size)


        if 'restrict_radius_pixels' in detconf:
            logger.debug('    restricting to radius %.2f '
                         'pixels' % detconf['restrict_radius_pixels'])
            # restrict to a central region
            ccen=(array(obs.image.shape)-1.0)/2.0
            rad=sqrt(
                (ccen[0]-objs['y'])**2
                +
                (ccen[1]-objs['x'])**2
            )
            w,=where(rad < detconf['restrict_radius_pixels'])
            logger.debug('    found %d within central region' % w.size)
            if w.size == 0:
                raise TryAgainError("no objects found in central region")
            objs = objs[w]

        logger.debug('    y: %s' % str(objs['y']))
        logger.debug('    x: %s' % str(objs['x']))
        if ('deblending' in detconf
                and detconf['deblending']['deblend']):
            #and objs.size > 1):
            if detconf['deblending']['type']=='full':
                list_of_obslist=self._run_mof(objs, obs)
            else:
                list_of_obslist=self._run_mof_stamps(objs, obs)
        else:
            list_of_obslist=self._extract_simple_stamps(objs, obs)

        if self['show']:
            self._show_images([o[0].image for o in list_of_obslist])

        return list_of_obslist

    def _run_mof_stamps(self, objs, obs):
        """
        fit all objects using MOF and get corrected observations

        Not yet doing FoF grouping
        """
        import mof

        # we need obs to have jacobian with center at 0,0 for the
        # full MOF
        debconf=self['detection']['deblending']
        mofconf=debconf['mof']

        # set the PSF

        pconf=mofconf['psf_pars']
        assert pconf['model']=='coellip2'
        jac=obs.jacobian
        Tguess=4.0*jac.scale**2
        ngauss=2
        runner=ngmix.bootstrap.PSFRunnerCoellip(
            obs.psf, 
            Tguess,
            ngauss,
            pconf['lm_pars'],
        )
        runner.go(ntry=pconf['ntry'])
        psf_fitter=runner.get_fitter()
        pres=psf_fitter.get_result()
        if pres['flags']!= 0:
            raise TryAgainError('psf fitting for MOF failed')

        psf_gmix=psf_fitter.get_gmix() 
        obs.psf.set_gmix(psf_gmix)

        # make the medsifier and MEDS object
        dlist=[
            {
                'image':obs.image,
                'weight':obs.weight,
                'wcs':obs.jacobian.get_galsim_wcs(),
            }
        ]
        mer=mof.stamps.MEDSifier(
            dlist,
            meds_config=debconf['meds_config'],
            sx_config=debconf['sx_config'],
        )
        mbm=mer.get_multiband_meds()
        list_of_mbobs=mbm.get_mbobs_list(weight_type=debconf['weight_type'])
        for mbobs in list_of_mbobs:
            mbobs[0][0].set_psf( obs.psf.copy() )

        prior=mof.moflib.get_mof_stamps_prior(
            list_of_mbobs,
            mofconf['model'],
            self.rng,
        )

        fitter=mof.MOFStamps(
            list_of_mbobs,
            mofconf['model'],
            prior=prior,
        )
        detband=0
        for i in xrange(mofconf['ntry']):
            guesses=mof.moflib.get_stamp_guesses(
                list_of_mbobs,
                detband,
                mofconf['model'],
                self.rng,
            )
            #ngmix.print_pars(guesses, front='    guess: ')
            fitter.go(guesses)
            res=fitter.get_result()
            if res['flags']==0:
                break

        if res['flags'] != 0:
            raise TryAgainError('mof failed')

        list_of_obslist = fitter.make_corrected_obs(band=0)
        #for obslist in list_of_obslist:
        #    assert obslist[0].has_psf_gmix()

        return list_of_obslist


    def _run_mof(self, objs, obs):
        """
        fit all objects using MOF and get corrected observations

        Not yet doing FoF grouping
        """
        import mof

        # we need obs to have jacobian with center at 0,0 for the
        # full MOF
        mofconf=self['detection']['deblending']['mof']

        jac=obs.jacobian
        jac.set_cen(row=0, col=0)
        logger.debug('jacobian: %s' % jac)
        obs.jacobian=jac

        pconf=mofconf['psf_pars']
        assert pconf['model']=='coellip2'
        Tguess=4.0*jac.scale**2
        ngauss=2
        runner=ngmix.bootstrap.PSFRunnerCoellip(
            obs.psf, 
            Tguess,
            ngauss,
            pconf['lm_pars'],
        )
        runner.go(ntry=pconf['ntry'])
        psf_fitter=runner.get_fitter()
        pres=psf_fitter.get_result()
        if pres['flags']!= 0:
            raise TryAgainError('psf fitting for MOF failed')

        psf_gmix=psf_fitter.get_gmix() 
        obs.psf.set_gmix(psf_gmix)

        nband=1
        prior=mof.moflib.get_mof_full_image_prior(
            objs,
            nband,
            jac,
            mofconf['model'],
            self.rng,
        )
        fitter=mof.MOF(
            obs,
            mofconf['model'],
            objs.size,prior=prior,
        )
        for i in xrange(mofconf['ntry']):
            guesses=mof.moflib.get_full_image_guesses(
                objs,
                nband,
                jac,
                mofconf['model'],
                self.rng,
                #Tguess=1.0*jac.scale**2,
            )
            #ngmix.print_pars(guesses, front='    guess: ')
            fitter.go(guesses)
            res=fitter.get_result()
            if res['flags']==0:
                break

        if res['flags'] != 0:
            raise TryAgainError('mof failed')

        list_of_obslist=[]
        for i in xrange(objs.size): 
            cobslist = fitter.make_corrected_obs(index=i, band=0)
            list_of_obslist.append(cobslist)

        return list_of_obslist

    def _extract_simple_stamps(self, objs, obs):
        # make list of obslist centered on these positions
        list_of_obslist=[]
        for i in xrange(objs.size):
            row=objs['y'][i]
            col=objs['x'][i]

            # this makes a copy
            jac=obs.jacobian

            jac.set_cen(row=row, col=col)

            newobs = obs.copy()
            newobs.jacobian=jac

            newobs = self._trim_image_from_detection(newobs)
            new_obslist=ngmix.ObsList()
            new_obslist.append(newobs)

            list_of_obslist.append( new_obslist )
        return list_of_obslist

    def _find_center_sep(self, obslist):
        from . import runsep

        assert len(obslist)==1
        obs=obslist[0]

        # this makes a copy
        jac=obs.jacobian

        cconf=self['center']

        objs, seg = runsep.find_objects(
            obs,
            segmentation_map=True,
        )

        if objs.size == 0:
            if self['show']:
                self._show_image_and_seg(obs.image, seg)
            raise TryAgainError('no objects found')

        logger.debug('    found %d objects' % objs.size)
        if objs.size > 1:
            logger.debug('        y:    %s' % (objs['y']))
            logger.debug('        x:    %s' % (objs['x']))
            logger.debug('        flux: %s' % (objs['flux']))

        if 'restrict_radius_pixels' in cconf:
            logger.debug('    restricting to radius %.2f '
                         'pixels' % cconf['restrict_radius_pixels'])
            # restrict to the central region
            ccen=(array(obs.image.shape)-1.0)/2.0
            rad=sqrt(
                (ccen[0]-objs['y'])**2
                +
                (ccen[1]-objs['x'])**2
            )
            w,=where(rad < cconf['restrict_radius_pixels'])
            logger.debug('    found %d within central region' % w.size)
            if w.size == 0:
                if self['show']:
                    logger.debug('no objects found in central region')
                    self._show_image_and_seg(obs.image, seg)
                raise TryAgainError("no objects found in central region")
            objs = objs[w]

        imax=objs['flux'].argmax()
        row=objs['y'][imax]
        col=objs['x'][imax]

        # update the jacobian center

        oldrow,oldcol=jac.get_cen()
        logger.debug('    old pos: %f %f' % (oldrow,oldcol))
        logger.debug('    new pos: %f %f' % (row,col))
        logger.debug('    shift: %f %f' % (row-oldrow,col-oldcol))

        crow,ccol=(array(obs.image.shape)-1.0)/2.0
        offset_from_ccen=(row-crow, col-ccol)
        logger.debug('    offset_from_ccen: %g %g ' % tuple(offset_from_ccen))

        if self['show']:
            self._show_image_and_seg(obs.image, seg)

        jac.set_cen(row=row, col=col)
        obs.jacobian=jac

    def _show_image(self, im):
        import images
        images.multiview(im)
        if 'q'==raw_input('hit a key (q to quit): '):
            stop

    def _show_images(self, imlist):
        import images
        images.view_mosaic(imlist)
        if 'q'==raw_input('hit a key (q to quit): '):
            stop



    def _show_image_and_seg(self, im, seg):
        import images
        images.view_mosaic([im, seg])
        if 'q'==raw_input('hit a key (q to quit): '):
            stop

    def _trim_image_from_detection(self, obs):

        dims = obs.image.shape

        jac = obs.jacobian

        # this is in arcsec
        cen=array(jac.get_cen())

        rowpix=int(round(cen[0]))
        colpix=int(round(cen[1]))

        # now get smallest distance to an edge
        radpix = self['detection']['stamp_radius']

        # note adding 1 for the slice
        row_start = rowpix-radpix
        row_end   = rowpix+radpix+1
        col_start = colpix-radpix
        col_end   = colpix+radpix+1
        logger.debug('    row start/end %s:%s' % (row_start, row_end))
        logger.debug('    col start/end %s:%s' % (col_start, col_end))

        im=obs.image
        subim = im[
            row_start:row_end,
            col_start:col_end,
        ]
        logger.debug('    subim dims: %s' % str(subim.shape))
        assert subim.shape[0]==subim.shape[1]

        newcen = cen - array([row_start, col_start])
        logger.debug('    new cen: %g %g' % tuple(newcen))
        new_ccen=(array(subim.shape)-1.0)/2.0
        logger.debug('    new ccen: %g %g' % tuple(new_ccen))

        offset_from_ccen=newcen-new_ccen
        logger.debug('    offset from ccen: %g %g ' % tuple(offset_from_ccen))

        jac.set_cen(row=newcen[0], col=newcen[1])

        wt00=obs.weight[0,0]
        newobs = ngmix.Observation(
            subim,
            weight=subim*0 + wt00,
            jacobian=jac,
            meta=obs.meta,
            psf=obs.psf.copy(),
        )
        return newobs

    def _trim_image_for_centering(self, obslist):
        assert len(obslist)==1
        obs=obslist[0]

        dims = obs.image.shape

        jac = obs.jacobian

        # this is in arcsec
        cen=array(jac.get_cen())

        rowpix=int(round(cen[0]))
        colpix=int(round(cen[1]))

        # now get smallest distance to an edge
        radpix = self['center'].get('trim_radius',None)
        if radpix is None:
            radpix = min(

                rowpix-0,
                (dims[0]-1)-rowpix,

                colpix-0,
                (dims[1]-1)-colpix,
            )

        # note adding 1 for the slice
        row_start = rowpix-radpix
        row_end   = rowpix+radpix+1
        col_start = colpix-radpix
        col_end   = colpix+radpix+1
        logger.debug('    row start/end %s:%s' % (row_start, row_end))
        logger.debug('    col start/end %s:%s' % (col_start, col_end))

        im=obs.image
        subim = im[
            row_start:row_end,
            col_start:col_end,
        ]
        logger.debug('    subim dims: %s' % str(subim.shape))
        assert subim.shape[0]==subim.shape[1]

        newcen = cen - array([row_start, col_start])
        logger.debug('    new cen: %g %g' % tuple(newcen))
        new_ccen=(array(subim.shape)-1.0)/2.0
        logger.debug('    new ccen: %g %g' % tuple(new_ccen))

        offset_from_ccen=newcen-new_ccen
        logger.debug('    offset from ccen: %g %g ' % tuple(offset_from_ccen))

        jac.set_cen(row=newcen[0], col=newcen[1])

        wt00=obs.weight[0,0]
        newobs = ngmix.Observation(
            subim,
            weight=subim*0 + wt00,
            jacobian=jac,
            meta=obs.meta,
            psf=obs.psf.copy(),
        )

        newobslist=ngmix.ObsList()
        newobslist.append(newobs)

        if self['show']:
            self._show_images([obs.image, newobs.image])

        return newobslist



    def _make_output(self, res):
        d=self._get_struct()
        return d

    def _setup(self, run_conf, **keys):
        """
        Check and set the configurations
        """

        self['s2n_min'] = self.get('s2n_min',-9999.0)
        self['deblend'] = self.get('deblend',False)

        self['nrand'] = self.get('nrand',1)
        if self['nrand'] is None:
            self['nrand']=1

        self.update(run_conf)

        self['use_logpars']=self.get('use_logpars',False)

        if 'fit_model' in self:
            if 'spergel' in self['fit_model']:
                self['npars']=7
            else:
                if 'galsim' in self['fitter']:
                    self['npars']=6
                else:
                    self['npars']=ngmix.gmix.get_model_npars(self['fit_model'])
        else:
            # default to 6
            self['npars']=6

        # this is "full" pars
        if 'psf_pars' in self:
            if 'am' in self['psf_pars']['model']:
                self['psf_npars']=6
            else:
                self['psf_npars']=\
                    6*ngmix.gmix._gmix_ngauss_dict[self['psf_pars']['model']]

        self['make_plots']=keys.get('make_plots',False)
        self['plot_base']=keys.get('plot_base',None)

        cenpars={
            'find_center':False,
            'find_center_metacal':False,
            'trim_image':False,
        }
        incenpars=self.get('center',{})
        cenpars.update(incenpars)
        self['center'] = cenpars
        if self['center']['find_center_metacal']:
            raise NotImplementedError('implement find_center_metacal')

    def _get_prior(self):
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

    def _get_dtype(self):
        """
        Set generic dtype
        """
        return []

    def _get_struct(self):
        """
        Make the output array
        """

        dt=self._get_dtype()
        data = numpy.zeros(1, dtype=dt)
        return data

class SimpleFitterBase(FitterBase):

    def _make_output(self, res):

        d = super(SimpleFitterBase,self)._make_output(res)

        n=self._get_namer()

        d[n('pars')] = res['pars']
        d[n('pars_cov')] = res['pars_cov']

        d['psf_pars'] = res['psf_pars']

        d[n('g')] = res['g']
        d[n('g_cov')] = res['g_cov']

        d[n('s2n_r')] = res['s2n_r']
        d['psf_T'] = res['psf_T']

        return d

    def _get_prior(self, prior_info=None):
        """
        Set all the priors
        """
        import ngmix
        from ngmix.joint_prior import PriorSimpleSep

        if 'priors' not in self:
            return None

        if 'masking' in self:
            self

        if prior_info is None:
            ppars=self['priors']
        else:
            ppars=prior_info

        gp = ppars['g']
        if gp['type']=='ba':
            g_prior = ngmix.priors.GPriorBA(gp['sigma'], rng=self.rng)
        elif gp['type']=='flat':
            g_prior = ngmix.priors.ZDisk2D(1.0, rng=self.rng)
        else:
            raise ValueError("implement other g prior")

        logger.debug("using input search prior")

        self.fracdev_prior=None
        if 'fracdev' in ppars:
            self.fracdev_prior = load_gmixnd(ppars['fracdev'], rng=self.rng)
            print("added fracdev prior:",self.fracdev_prior)

        T_prior = None
        if 'T' in ppars:
            Tp = ppars['T']

            if Tp['type']=="flat":
                T_prior=ngmix.priors.FlatPrior(*Tp['pars'], rng=self.rng)

            elif Tp['type']=="gmixnd":
                T_prior = load_gmixnd(Tp, rng=self.rng)

            elif Tp['type']=='normal':
                Tpars=Tp['pars']
                T_prior=ngmix.priors.Normal(
                    Tp['mean'],
                    Tp['sigma'],
                    rng=self.rng,
                )

            elif Tp['type']=='lognormal':

                if self['use_logpars']:
                    logger.debug("    converting T prior to log")
                    logT_mean, logT_sigma=ngmix.priors.lognorm_convert(Tp['mean'],
                                                                       Tp['sigma'])
                    T_prior = ngmix.priors.Normal(logT_mean, logT_sigma, rng=self.rng)

                else:
                    shift=Tp.get('shift',None)
                    T_prior = ngmix.priors.LogNormal(
                        Tp['mean'],
                        Tp['sigma'],
                        shift=shift,
                        rng=self.rng,
                    )

            elif Tp['type'] in ['TwoSidedErf',"two-sided-erf"]:
                T_prior_pars = Tp['pars']
                T_prior=ngmix.priors.TwoSidedErf(*T_prior_pars, rng=self.rng)
            else:
                raise ValueError("bad Tprior: '%s'" % Tp['type'])


        r50_prior = None
        if 'r50' in ppars:
            assert self['use_logpars']==False

            r50p = ppars['r50']

            if r50p['type']=="gmixnd":
                r50_prior = load_gmixnd(r50p, rng=self.rng)

            elif r50p['type']=='normal':
                T_prior=ngmix.priors.Normal(
                    r50p['mean'],
                    r50p['sigma'],
                    rng=self.rng,
                )


            elif r50p['type']=='lognormal':

                r50_prior = ngmix.priors.LogNormal(
                    r50p['mean'],
                    r50p['sigma'],
                    shift=r50p.get('shift',None),
                    rng=self.rng,
                )

            elif r50p['type'] in ['TwoSidedErf',"two-sided-erf"]:
                r50_prior=ngmix.priors.TwoSidedErf(*r50p['pars'], rng=self.rng)

            elif r50p['type']=="flat":
                r50_prior=ngmix.priors.FlatPrior(*r50p['pars'], rng=self.rng)
            else:
                raise ValueError("bad r50prior: '%s'" % r50p['type'])

        nu_prior=None
        if 'nu' in ppars:

            nup = ppars['nu']

            if nup['type']=="gmixnd":
                nu_prior = load_gmixnd(nup, rng=self.rng)

            elif nup['type']=='lognormal':
                nu_prior = ngmix.priors.LogNormal(
                    nup['mean'],
                    nup['sigma'],
                    shift=nup.get('shift',None),
                    rng=self.rng,
                )

            elif nup['type']=='normal':
                nu_prior = ngmix.priors.Normal(
                    nup['mean'],
                    nup['sigma'],
                    rng=self.rng,
                )


            elif nup['type'] in ['TwoSidedErf',"two-sided-erf"]:
                nu_prior=ngmix.priors.TwoSidedErf(*nup['pars'], rng=self.rng)

            elif nup['type']=="flat":
                nu_prior=ngmix.priors.FlatPrior(*nup['pars'], rng=self.rng)
            else:
                raise ValueError("bad nuprior: '%s'" % nup['type'])

        n_prior=None
        if 'n' in ppars:

            np = ppars['n']

            if np['type']=="gmixnd":
                n_prior = load_gmixnd(np, rng=self.rng)

            elif np['type']=='lognormal':
                n_prior = ngmix.priors.LogNormal(
                    np['mean'],
                    np['sigma'],
                    shift=np.get('shift',None),
                    rng=self.rng,
                )

            elif np['type']=='normal':
                n_prior = ngmix.priors.Normal(
                    np['mean'],
                    np['sigma'],
                    rng=self.rng,
                )


            elif np['type'] in ['TwoSidedErf',"two-sided-erf"]:
                n_prior=ngmix.priors.TwoSidedErf(*np['pars'], rng=self.rng)

            elif np['type']=="flat":
                n_prior=ngmix.priors.FlatPrior(*np['pars'], rng=self.rng)
            else:
                raise ValueError("bad nprior: '%s'" % np['type'])



        cp=ppars['counts']
        if cp['type']=="truth":
            logger.debug("using true counts pdf for prior")

            if self['use_logpars']:

                logger.debug("    converting to log")
                if self.sim['simulator']=="galsim":
                    fluxspec=self.sim['obj_model']['flux']
                    counts       = fluxspec['mean']
                    counts_sigma = fluxspec['sigma']
                else:
                    counts       = self.sim['obj_counts_mean']
                    counts_sigma = self.sim['obj_counts_sigma_frac']*counts

                logc_mean, logc_sigma=ngmix.priors.lognorm_convert(counts,
                                                                   counts_sigma)
                counts_prior = ngmix.priors.Normal(logc_mean, logc_sigma, rng=self.rng)

            else:
                counts_prior = self.sim.counts_pdf

        elif cp['type']=="gmixnd":

            counts_prior = load_gmixnd(cp, rng=self.rng)

        elif cp['type']=='lognormal':

            if self['use_logpars']:
                logger.debug("    converting counts prior to log")
                logcounts_mean, logcounts_sigma=ngmix.priors.lognorm_convert(cp['mean'],
                                                                   cp['sigma'])
                counts_prior = ngmix.priors.Normal(
                    logcounts_mean,
                    logcounts_sigma,
                    rng=self.rng,
                )

            else:
                counts_prior = ngmix.priors.LogNormal(
                    cp['mean'],
                    cp['sigma'],
                    shift=cp.get('shift',None),
                    rng=self.rng,
                )




        elif cp['type']=='normal':
            counts_prior=ngmix.priors.Normal(
                cp['mean'],
                cp['sigma'],
                rng=self.rng,
            )

        elif cp['type'] in ['TwoSidedErf',"two-sided-erf"]:
            counts_prior=ngmix.priors.TwoSidedErf(*cp['pars'], rng=self.rng)

        elif cp['type']=="flat":
            counts_prior=ngmix.priors.FlatPrior(*cp['pars'], rng=self.rng)

        else:
            raise ValueError("bad counts prior: '%s'" % cp['type'])

        cp=ppars['cen']
        if cp['type']=="truth":
            cen_prior=self.sim.cen_pdf
        elif cp['type'] == "normal2d":
            fit_cen_sigma=cp['sigma']
            cen_prior=ngmix.priors.CenPrior(
                0.0,
                0.0,
                fit_cen_sigma,
                fit_cen_sigma,
                rng=self.rng,
            )
        else:
            raise ValueError("bad cen prior: '%s'" % cp['type'])

        if n_prior is not None:
            prior=ngmix.joint_prior.PriorSersicSep(
                cen_prior,
                g_prior,
                T_prior,
                n_prior,
                counts_prior,
            )
        else:

            if nu_prior is not None:
                prior = ngmix.joint_prior.PriorSpergelSep(
                    cen_prior,
                    g_prior,
                    r50_prior,
                    nu_prior,
                    counts_prior,
                )
            elif r50_prior is not None:
                prior = PriorSimpleSep(cen_prior,
                                       g_prior,
                                       r50_prior,
                                       counts_prior)

            else:
                prior = PriorSimpleSep(cen_prior,
                                       g_prior,
                                       T_prior,
                                       counts_prior)
        return prior


    def _get_namer(self):
        if 'fit_model' in self:
            front=self['fit_model']
        else:
            front=None

        return util.Namer(front=front)

    def _get_dtype(self):
        """
        get the dtype for the output struct
        """
        npars=self['npars']

        n=self._get_namer()

        dt=super(SimpleFitterBase,self)._get_dtype()
        dt += [
            (n('pars'),     'f8',npars),
            (n('pars_cov'), 'f8',(npars,npars)),
            (n('g'),        'f8',2),
            (n('g_cov'),    'f8',(2,2)),
            (n('s2n_r'),    'f8'),
        ]

        if 'psf_npars' in self:
            psf_npars=self['psf_npars']
            dt += [
                ('psf_pars','f8',psf_npars),
                ('psf_T','f8')
            ]

        return dt

class MaxFitter(SimpleFitterBase):

    def _dofit(self, obslist):
        """
        Fit according to the requested method
        """

        if 'mof' in self:
            obslist = self._do_mof_fit(obslist)

        use_round_T=self['use_round_T']
        if self['fit_model'] == 'cm':
            boot=ngmix.CompositeBootstrapper(
                obslist,
                fracdev_prior = self.fracdev_prior,
                verbose=False,
            )
        else:
            boot=ngmix.Bootstrapper(obslist,
                                    use_logpars=self['use_logpars'],
                                    use_round_T=use_round_T,
                                    verbose=False)

        mconf=self['max_pars']
        covconf=mconf['cov']

        Tguess=self._get_psf_T_guess(obslist)

        ppars=self['psf_pars']
        try:
            boot.fit_psfs(
                ppars['model'],
                Tguess,
                ntry=mconf['ntry'],
            )
        except BootPSFFailure:
            raise TryAgainError("failed to fit psf")


        try:

            boot.fit_max(self['fit_model'],
                         mconf['pars'],
                         prior=self.prior,
                         ntry=mconf['ntry'])

            boot.set_round_s2n()
            rres=boot.get_round_result()
            res=boot.get_max_fitter().get_result()

            res['psf_pars'] = boot.mb_obs_list[0][0].psf.gmix.get_full_pars()

            res['s2n_r'] = rres['s2n_r']
            res['T_r'] = rres['T_r']

            res['psf_T'] = obslist[0].psf.gmix.get_T()
            res['psf_T_r'] = rres['psf_T_r']

            pres=boot.get_psf_flux_result()
            res['psf_flux'] = pres['psf_flux']
            res['psf_flux_err'] = pres['psf_flux_err']

            if covconf['replace_cov']:
                boot.try_replace_cov(covconf['cov_pars'])

            if self['use_round_T']:
                res['T_s2n_r'] = boot.get_max_fitter().get_T_s2n()

        except BootGalFailure:
            raise TryAgainError("failed to fit galaxy")

        fitter=boot.get_max_fitter() 

        return fitter

    def _do_mof_fit(self, allobs):
        import minimof

        try:
            mm=minimof.MiniMOF(
                self['mof'],
                allobs,
                self.rng,
            )
            mm.go()
            mm_res = mm.get_result()
            if not mm_res['converged']:
                raise TryAgainError("MOF did not converge")

            if False:
                mm.show_residuals()
                if raw_input("q to quit: ")=='q':
                    stop
        except BootPSFFailure as err:
            raise TryAgainError("MOF psf failure: '%s'" % str(err))
        except BootGalFailure as err:
            raise TryAgainError("MOF gal failure: '%s'" % str(err))

        # assume first is the central
        corr_obslist = mm.get_corrected_obs(0)

        return corr_obslist


    def _get_psf_T_guess(self, obslist):
        if isinstance(obslist,ngmix.observation.ObsList):
            j=obslist[0].jacobian
        else:
            j=obslist.jacobian

        scale=j.get_scale()
        Tguess = 4.0*scale**2
        return Tguess
 
    def _print_res(self,res):
        """
        print some stats
        """

        if 'nfev' in res:
            mess="    s2n: %.1f  ntry: %d  nfev: %d"
            mess = mess % (res['s2n'],res['ntry'],res['nfev'])
            logger.debug(mess)

        log_pars(res['pars'],      front='        pars: ')
        log_pars(res['pars_err'],  front='        perr: ')

    def _make_plots(self, fitter, key):
        """
        Write a plot file of the trials
        """
        import biggles

        biggles.configure('default','fontsize_min',1.0)

        width,height=1100,1100

        #pdict=fitter.make_plots(do_residual=True, title=self.fit_model)

        resid_pname=self['plot_base']+'-%06d-%s-gal-resid.png' % (self.igal,key)
        logger.debug(resid_pname)
        rplt=fitter.plot_residuals()
        rplt[0][0].write_img(width,height,resid_pname)


    def _get_dtype(self):
        """
        get the dtype for the output struct
        """
        npars=self['npars']

        dt=super(MaxFitter,self)._get_dtype()
        dt += [
            ('s2n_r','f8'),
            ('T_r','f8'),
            ('psf_T_r','f8'),
            ('nfev','i4'),
            ('ntry','i4'),
            ('psf_flux','f8'),
            ('psf_flux_err','f8'),
        ]

        if self['use_round_T']:
            dt += [('T_s2n_r','f8')]

        return dt


    def _make_output(self, res):

        d = super(MaxFitter,self)._make_output(res)

        if 'nfev' in res:
            d['s2n_r'] = res['s2n_r']
            d['T_r'] = res['T_r']
            d['psf_T_r'] = res['psf_T_r']
            d['nfev'] = res['nfev']
            # set outside of fitter
            d['ntry'] = res['ntry']

            if 'psf_flux' in res:
                d['psf_flux'] = res['psf_flux']
                d['psf_flux_err'] = res['psf_flux_err']

            if self['use_round_T']:
                d['T_s2n_r'] = res['T_s2n_r']

        return d

class MaxMetacalFitter(MaxFitter):
    """
    metacal with a maximum likelihood fit
    """
    def __init__(self, *args, **kw):
        super(MaxMetacalFitter,self).__init__(*args, **kw)

        self['min_s2n'] = self.get('min_s2n',0.0)

    def _dofit(self, obslist):
        """
        Fit according to the requested method
        """

        if 'mof' in self:
            obslist = self._do_mof_fit(obslist)

        mdict = self._do_fits(obslist)
        res=mdict['res']
        fitter=mdict['fitter']

        self.fitter=fitter

        return fitter

        
    def _get_bootstrapper(self, obs):
        from ngmix.bootstrap import MaxMetacalBootstrapper
        metacal_pars=self['metacal_pars']
        if metacal_pars.get('shear_psf_inv',False):
            boot=MaxMetacalBootstrapperShearPSFInv(
                obs,
                verbose=False,
            )

        else:

            boot=MaxMetacalBootstrapper(obs,
                                     use_logpars=self['use_logpars'],
                                     verbose=False)
        return boot


    def _do_fits(self, obs):
        """
        the basic fitter for this class
        """

        if 'masking' in self:
            replace_fitter=self._do_fits_for_replacement(obs)

        if 'center' in self:
            if self['center']['find_center']:
                logger.debug("    finding center")
                self._find_center(obs)

            if self['center']['trim_image']:
                obs=self._trim_image_for_centering(obs)

        boot=self._get_bootstrapper(obs)

        mconf=self['max_pars']
        covconf=mconf['cov']

        ppars=self['psf_pars']

        try:
            Tguess=self._get_psf_T_guess(obs)
            # redo psf in case we did replacement fit above
            boot.fit_psfs(
                ppars['model'],
                Tguess,
                ntry=mconf['ntry'],
                fit_pars=mconf['pars']['lm_pars'],
                skip_already_done=False,
            )

        except BootPSFFailure:
            raise TryAgainError("failed to fit psf")


        try:
            boot.fit_max(self['fit_model'],
                         mconf['pars'],
                         prior=self.prior,
                         ntry=mconf['ntry'])

            boot.set_round_s2n()
            rres=boot.get_round_result()
            res=boot.get_max_fitter().get_result()

            res['psf_pars'] = boot.mb_obs_list[0][0].psf.gmix.get_full_pars()

            res['s2n_r'] = rres['s2n_r']

            if res['s2n_r'] < self['min_s2n']:
                raise TryAgainError("s2n %g < %g" % (res['s2n_r'],self['min_s2n']))

            res['T_r'] = rres['T_r']

            #res['psf_T'] = obs.psf.gmix.get_T()
            res['psf_T_r'] = rres['psf_T_r']
            res['psf_T'] = rres['psf_T_r']

            if covconf['replace_cov']:
                boot.try_replace_cov(covconf['cov_pars'])

            if 'masking' in self:
                boot.replace_masked_pixels(fitter=replace_fitter)

            self._do_metacal(boot)

        except BootPSFFailure:
            raise TryAgainError("failed to fit metacal psfs")
        except BootGalFailure:
            raise TryAgainError("failed to fit galaxy")

        fitter=boot.get_max_fitter()
        res=fitter.get_result()

        mres = boot.get_metacal_result()
        res.update(mres)

        return {'fitter':fitter,
                'boot':boot,
                'res':res}

    def _do_metacal(self, boot):

        ppars=self['psf_pars']
        mconf=self['max_pars']
        Tguess=self.sim.get('psf_T',4.0)
        psf_fit_pars = ppars.get('fit_pars',None)

        boot.fit_metacal(
            ppars['model'],
            self['fit_model'],
            mconf['pars'],
            Tguess,
            psf_fit_pars=psf_fit_pars,
            prior=self.prior,
            ntry=mconf['ntry'],
            metacal_pars=self['metacal_pars'],
        )

    def _find_center(self, obslist):
        if self['center']['method']=='sep':
            self._find_center_sep(obslist)
        else:
            raise ValueError("only sep used for finding center for now")


    def _get_dtype(self):
        """
        get the dtype for the output struct
        """
        dt=super(MaxMetacalFitter,self)._get_dtype()

        npars=self['npars']

        mpars=self['metacal_pars']
        types=mpars.get('types',ngmix.metacal.METACAL_TYPES)

        for t in ngmix.metacal.METACAL_REQUIRED_TYPES:
            if t not in types:
                types.append(t)

        for type in types:

            if type=='noshear':
                back=''
            else:
                back='_%s' % type

            dt += [
                ('mcal_g%s' % back,'f8',2),
                ('mcal_g_cov%s' % back,'f8',(2,2)),
                ('mcal_pars%s' % back,'f8',npars),
            ]

            if type=='noshear':
                dt += [
                    ('mcal_pars_cov','f8',(npars,npars)),
                    ('mcal_psfrec_g','f8',2),
                    ('mcal_psfrec_T','f8'),
                ]

            dt += [
                ('mcal_s2n_r%s' % back,'f8'),
                ('mcal_T_r%s' % back,'f8'),
            ]

        return dt


    def _make_output(self, res):
        """
        copy parameters specific to this class
        """
        d = super(MaxMetacalFitter,self)._make_output(res)


        for type in ngmix.metacal.METACAL_TYPES:

            # sometimes we don't calculate all
            if type not in res:
                continue

            tres=res[type]
            if type=='noshear':
                back=''
            else:
                back='_%s' % type

            d['mcal_pars%s' % back] = tres['pars']
            d['mcal_g%s' % back] = tres['g']
            d['mcal_g_cov%s' % back] = tres['g_cov']
            d['mcal_s2n_r%s' % back] = tres['s2n_r']
            d['mcal_T_r%s' % back] = tres['T_r']

            if type=='noshear':
                for p in ['pars_cov','psfrec_g','psfrec_T']:
                    name='mcal_%s' % p

                    if p=='psfrec_g':
                        d[name] = tres['gpsf']
                    elif p=='psfrec_T':
                        d[name] = tres['Tpsf']
                    else:
                        d[name] = tres[p]
        return d

    def _print_res(self,resfull):
        """
        print some stats
        """

        super(MaxMetacalFitter,self)._print_res(resfull)

        res=resfull['noshear']
        logger.debug("    mcal s2n_r: %g" % res['s2n_r'])
        log_pars(res['pars'],       front='    mcal pars: ')

class MaxMetacalBootstrapperShearPSFInv(ngmix.bootstrap.MaxMetacalBootstrapper):
    def _get_all_metacal(self, metacal_pars, **kw):
        metacal_pars=self._extract_metacal_pars(metacal_pars)
        return get_all_metacal_shearpsfinv(self.mb_obs_list, **metacal_pars)


class MetacalShearPSFInv(ngmix.metacal.Metacal):
    def get_target_psf(self, shear, type, get_nopix=False):
        """
        get galsim object for the dilated, possibly sheared, psf

        parameters
        ----------
        shear: ngmix.Shape
            The applied shear
        type: string
            Type of psf target.  For type='gal_shear', the psf is just dilated to
            deal with noise amplification.  For type='psf_shear' the psf is also
            sheared for calculating Rpsf

        returns
        -------
        galsim object
        """

        ngmix.metacal._check_shape(shear)
        
        denom = galsim.Deconvolve(self.psf_int_nopix).shear(
            g1=shear.g1,
            g2=shear.g2,
        )
        psf_grown_nopix = galsim.Deconvolve(denom)
        psf_grown = galsim.Convolve(psf_grown_nopix,self.pixel)

        # this should carry over the wcs
        psf_grown_image = self.psf_image.copy()

        try:
            psf_grown.drawImage(
                image=psf_grown_image,
                method='no_pixel' # pixel is in the psf
            )
            import images
            plt=images.multiview(psf_grown_image.array)
            plt.write('shc.png',dpi=150)
            stop

            if get_nopix:
                psf_grown_nopix_image = self.psf_image.copy()
                psf_grown_nopix.drawImage(
                    image=psf_grown_nopix_image,
                    method='no_pixel' # pixel is in the psf
                )

                return psf_grown_image, psf_grown_nopix_image, psf_grown
            else:
                return psf_grown_image, psf_grown


        except RuntimeError as err:
            # argh, galsim uses generic exceptions
            raise GMixRangeError("galsim error: '%s'" % str(err))


def get_all_metacal_shearpsfinv(obs, step=0.01, **kw):
    """
    internal routine

    get all metacal
    """
    if isinstance(obs, Observation):

        m=MetacalShearPSFInv(obs, **kw)
        odict=m.get_all(step, **kw)
    elif isinstance(obs, ngmix.MultiBandObsList):
        odict=make_metacal_mb_obs_list_dict_shearpsfinv(obs, step, **kw)
    elif isinstance(obs, ngmix.ObsList):
        odict=make_metacal_obs_list_dict_shearpsfinv(obs, step, **kw)
    else:
        raise ValueError("obs must be Observation, ObsList, "
                         "or MultiBandObsList")

    return odict

def make_metacal_mb_obs_list_dict_shearpsfinv(mb_obs_list, step, **kw):

    new_dict=None
    for obs_list in mb_obs_list:
        odict = make_metacal_obs_list_dict_shearpsfinv(obs_list, step, **kw)

        if new_dict is None:
            new_dict=ngmix.metacal._init_mb_obs_list_dict(odict.keys())

        for key in odict:
            new_dict[key].append(odict[key])

    return new_dict


def make_metacal_obs_list_dict_shearpsfinv(obs_list, step, **kw):
    odict = None
    first=True
    for obs in obs_list:

        todict=get_all_metacal_shearpsfinv(obs, step=step, **kw)

        if odict is None:
            odict=ngmix.metacal._init_obs_list_dict(todict.keys())

        for key in odict:
            odict[key].append( todict[key] )

    return odict


class GalsimFitter(SimpleFitterBase):

    def _get_guesser(self, obs):
        r50guess = 1.0
        flux_guess = obs.image.sum()

        guesser=ngmix.guessers.R50FluxGuesser(
            r50guess,
            flux_guess,
            prior=self.prior,
        )

        #guesser=ngmix.guessers.PriorGuesser(self.prior)   

        return guesser

    def _get_runner(self, obs):

        guesser=self._get_guesser(obs)

        mconf=self['max_pars']
        runner=ngmix.galsimfit.GalsimRunner(
            obs,
            self['fit_model'],
            guesser,
            lm_pars=mconf['lm_pars'],
            prior=self.prior,
        )
        return runner

    def _fit_am(self, obs, ntry=4):
        rng=self.rng

        am=ngmix.admom.Admom(obs, maxiter=1000)
        
        scale=obs.jacobian.get_scale()
        Tguess=4.0*scale

        grange=0.02

        for i in xrange(ntry):
            pars=[
                rng.uniform(low=-0.1*scale, high=0.1*scale),
                rng.uniform(low=-0.1*scale, high=0.1*scale),
                rng.uniform(low=-grange,high=grange),
                rng.uniform(low=-grange,high=grange),
                Tguess*(1.0 + rng.uniform(low=-0.1,high=0.1)),
                1.0,
            ]
            guess=ngmix.GMixModel(pars, "gauss")

            am.go(guess)
            res=am.get_result()
            if res['flags']==0:
                break

        if res['flags'] != 0:
            raise TryAgainError("failed to fit psf: %s" % pres['flagstr'])

        pgm=am.get_gmix()
        g1,g2,T=pgm.get_g1g2T()

        return g1,g2,T

    def _dofit(self, obslist):
        """
        Fit according to the requested method
        """

        assert len(obslist)==1
        obs = obslist[0]

        mconf=self['max_pars']

        try:
            runner=self._get_runner(obs)

            runner.go(ntry=mconf['ntry'])

            fitter=runner.get_fitter() 
            res=fitter.get_result()
            if res['flags'] != 0:
                raise TryAgainError("failed to fit galaxy")

            g1,g2,T=self._fit_am(obs.psf)
            res['gpsf'] = numpy.array([g1,g2])
            res['Tpsf'] = T

        except GMixRangeError as err:
            raise TryAgainError("failed to fit galaxy: %s" % str(err))

        return fitter

    def _print_res(self,res):
        """
        print some stats
        """

        if 'nfev' in res:
            mess="    s2n_r: %.1f  ntry: %d  nfev: %d"
            mess = mess % (res['s2n_r'],res['ntry'],res['nfev'])
            logger.debug(mess)

        log_pars(res['pars'],      front='        pars: ')
        log_pars(res['pars_err'],  front='        perr: ')

        #if res['pars_true'][0] is not None:
        #    log_pars(res['pars_true'], front='        true: ')


    def _get_dtype(self):
        """
        get the dtype for the output struct
        """
        npars=self['npars']

        dt=super(GalsimFitter,self)._get_dtype()
        dt += [
            ('s2n_r','f8'),
            ('nfev','i4'),
            ('ntry','i4')
        ]

        return dt

    def _make_output(self, res, i):

        # note super here
        d = super(SimpleFitterBase,self)._make_output(res)

        n=self._get_namer()
        d[n('pars')][i,:] = res['pars']
        d[n('pars_cov')][i,:,:] = res['pars_cov']

        if 'psf_pars' in res:
            d['psf_pars'][i,:] = res['psf_pars']

        d[n('g')][i,:] = res['g']
        d[n('g_cov')][i,:,:] = res['g_cov']

        if 'psf_T' in res:
            d['psf_T'] = res['psf_T']

        if 's2n_r' in res:
            d[n('s2n_r')] = res['s2n_r']

        if 'nfev' in res:
            d['nfev'] = res['nfev']
            # set outside of fitter
            d['ntry'] = res['ntry']

        return d

class SpergelFitter(GalsimFitter):

    def _get_guesser(self, obs):
        r50guess_pixels = 2.0
        scale=obs.jacobian.get_scale()
        r50guess = r50guess_pixels*scale
        flux_guess = obs.image.sum()

        # 0.5 is equivalent to sersic n=1, exponential
        nuguess=numpy.random.uniform(low=0.0,high=1.0)
        #nuguess = self.prior.nu_prior.sample()

        if True==self.get('guess_prior',False):
            guesser=ngmix.guessers.PriorGuesser(self.prior)
        else:
            guesser=ngmix.guessers.R50NuFluxGuesser(
                r50guess,
                nuguess,
                flux_guess,
                prior=self.prior,
            )

        return guesser






class SpergelMetacalFitter(SpergelFitter):
    def _setup(self, *args, **kw):
        super(SpergelMetacalFitter,self)._setup(*args, **kw)

        self['metacal_pars'] = self.get('metacal_pars',{})

        mpars=self['metacal_pars']
        self.metacal_types=mpars.get('types',ngmix.metacal.METACAL_TYPES)

        for t in ngmix.metacal.METACAL_REQUIRED_TYPES:
            if t not in self.metacal_types:
                self.metacal_types.append(t)

    def _dofit(self, imdict):
        """
        Fit according to the requested method
        """


        obs=imdict['obs']
        mconf=self['max_pars']

        if mconf['cov']['replace_cov']:
            replace_cov=True
            cov_pars=mconf['cov']['cov_pars']
        else:
            replace_cov=False

        res={}
        try:

            mc=ngmix.metacal.Metacal(obs)
            mcpars=self['metacal_pars']


            odict=ngmix.metacal.get_all_metacal(
                obs,
                rng=self.rng,
                **mcpars
            )

            for type in odict:
                mobs=odict[type]
                # for now only a single obs, since below we assume so to fit
                # the gpsf
                assert isinstance(mobs,Observation)

                runner=self._get_runner(mobs)

                runner.go(ntry=mconf['ntry'])

                fitter=runner.get_fitter() 
                
                tres=fitter.get_result()
                if tres['flags'] != 0:
                    raise TryAgainError("failed to fit a metacal obs")

                if replace_cov:
                    fitter.calc_cov(cov_pars['h'],cov_pars['m'])

                    if tres['flags'] != 0:
                        logger.debug("        cov replacement failed")
                        tres['flags']=0


                if type=='noshear':
                    g1,g2,T=self._fit_am(mobs.psf_nopix)
                    tres['gpsf'] = numpy.array([g1,g2])
                    tres['Tpsf'] = T

                res[type] = tres


        except GMixRangeError as err:
            raise TryAgainError("failed to fit galaxy: %s" % str(err))

        # just return the last fitter used
        res['flags']=0
        return res


    def _get_dtype(self):
        """
        get the dtype for the output struct
        """
        npars=self['npars']

        # super of super
        dt=super(SpergelFitter,self)._get_dtype()

        for type in self.metacal_types:

            if type=='noshear':
                back=''
            else:
                back='_%s' % type

            dt += [
                ('mcal_g%s' % back,'f8',2),
                ('mcal_g_cov%s' % back,'f8',(2,2)),
                ('mcal_pars%s' % back,'f8',npars),
            ]

            if type=='noshear':
                dt += [
                    ('mcal_pars_cov','f8',(npars,npars)),
                    ('mcal_gpsf','f8',2),
                    ('mcal_Tpsf','f8'),
                ]

            dt += [
                ('mcal_s2n_r%s' % back,'f8'),

                ('mcal_r50%s' % back,'f8'),
                ('mcal_r50_s2n%s' % back,'f8'),
                ('mcal_flux%s' % back,'f8'),
                ('mcal_flux_s2n%s' % back,'f8'),
            ]

        return dt

    def _make_output(self, res, i):
        """
        copy parameters specific to this class
        """

        # note copying super of our super, since
        # we didn't do a regular fit
        d = super(SimpleFitterBase,self)._make_output(res)

        for type in self.metacal_types:

            tres=res[type]
            if type=='noshear':
                back=''
            else:
                back='_%s' % type

            d['mcal_pars%s' % back] = tres['pars']
            d['mcal_g%s' % back] = tres['g']
            d['mcal_g_cov%s' % back] = tres['g_cov']

            d['mcal_s2n_r%s' % back] = tres['s2n_r']

            r50 = tres['pars'][4]
            r50_s2n = r50/sqrt(tres['pars_cov'][4,4])
            d['mcal_r50%s' % back] = r50
            d['mcal_r50_s2n%s' % back] = r50_s2n

            flux     = tres['pars'][5]
            flux_s2n = flux/sqrt(tres['pars_cov'][5,5])
            d['mcal_flux%s' % back] = flux
            d['mcal_flux_s2n%s' % back] = flux_s2n

            if type=='noshear':
                for p in ['pars_cov','gpsf','Tpsf']:

                    if p in tres:
                        name='mcal_%s' % p
                        d[name] = tres[p]
        return d

    def _print_res(self,res):
        """
        print some stats
        """

        subres=res['noshear']

        mess="    mcal s2n_r: %.1f nfev: %d"
        logger.debug(mess % (subres['s2n_r'],subres['nfev']))


        log_pars(subres['pars'],      front='        pars: ')
        log_pars(subres['pars_err'],  front='        perr: ')

        log_pars(res['pars_true'], front='        true: ')


class GalsimMetacalFitter(SpergelMetacalFitter):
    # argh, need to fix this, the code is duplicated
    def _get_guesser(self, obs):

        r50guess_pixels = 2.0
        scale=obs.jacobian.get_scale()
        r50guess = r50guess_pixels*scale
        flux_guess = obs.image.sum()

        guesser=ngmix.guessers.R50FluxGuesser(
            r50guess,
            flux_guess,
            prior=self.prior,
        )

        return guesser

    def _get_runner(self, obs):

        guesser=self._get_guesser(obs)

        mconf=self['max_pars']
        runner=ngmix.galsimfit.GalsimRunner(
            obs,
            self['fit_model'],
            mconf['lm_pars'],
            guesser,
            prior=self.prior,
        )
        return runner


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


