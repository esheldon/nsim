"""
moment based code
"""
from __future__ import print_function
try:
    xrange
except:
    xrange=range
    raw_input=input

import logging
import numpy
from numpy import array, zeros, ones, log, log10, exp, sqrt, diag
from numpy import where, isfinite

import ngmix
from .util import log_pars
from ngmix.shape import Shape
from ngmix.gexceptions import BootPSFFailure, BootGalFailure

from . import util
from .util import TryAgainError

from .fitters import SimpleFitterBase

import galsim

logger = logging.getLogger(__name__)

DEFTYPES=[
    'noshear',
    '1p','1m','2p','2m',
]



class MetacalMomentsAM(SimpleFitterBase):
    """
    TODO:

        - flux no longer from am

    """
    def _setup(self, *args, **kw):
        super(MetacalMomentsAM,self)._setup(*args, **kw)

        self['metacal_pars'] = self.get('metacal_pars',{})

        mpars=self['metacal_pars']

        if 'types' not in mpars:
            deftypes=[
                'noshear',
                '1p','1m','2p','2m',
            ]
            sym=mpars.get('symmetrize_psf',False)
            if not sym and 'psf' not in mpars:
                deftypes += [
                    '1p_psf','1m_psf',
                    '2p_psf','2m_psf',
                ]

                if 'shear_pixelized_psf' in mpars:
                    assert mpars['shear_pixelized_psf']==True
                elif 'prepix' in mpars:
                    assert mpars['prepix']==True
                else:
                    raise ValueError("if not symmetrizing with am, must "
                                     "set shear_pixelized_psf "
                                     "or prepix")
            mpars['types'] = deftypes

        self.metacal_types=mpars['types']
        logger.debug("doing types: %s" % self.metacal_types)

        self._set_mompars()

        self['min_s2n'] = self.get('min_s2n',0.0)

    def _set_mompars(self):
        self.ampars = self.get('admom_pars',{})
        self.psf_ampars = {}
        self.psf_ampars.update(self.ampars)
        self.psf_ampars['fixcen']=False
        self.psf_ampars['round']=False
        self.psf_ampars['use_canonical_center']=False

    def _dofit(self, obslist):


        if 'fit_model' in self:
            # measuring fluxes.  We will use this model for template
            # flux fitting in the metacal step
            flux_fitter = self._fit_flux_model(obs)
            flux_result=flux_fitter.get_result()
            if flux_result['flags'] != 0:
                raise TryAgainError("        flux fitting failed")

            gpars = flux_result['pars'].copy()
            gpars[0:0+2] = 0.0
            self.model_galsim_obj = flux_fitter.make_model(gpars)

        psfres = self._measure_psfs(obslist)

        obsdict=self._get_metacal(obslist)

        res=self._do_metacal(obsdict)
        res['psf'] = psfres


        if 'fit_model' in self:
            res['model_result'] = flux_result

        res['flags']=0
        return res

    def _measure_psfs(self, obslist):
        """
        measure the mean psf size and shape
        """
        Tsum=0.0
        g1sum=0.0
        g2sum=0.0

        for obs in obslist:
            psfres, fitter = self._measure_psf_moments(obs.psf)
            Tsum  += psfres['T']
            g1sum += psfres['g'][0]
            g2sum += psfres['g'][1]

            if fitter is not None:
                gmix=fitter.get_gmix()
                obs.psf.set_gmix(gmix)

        n=len(obslist)
        return {
            'T':Tsum/n,
            'g': [g1sum/n, g2sum/n],
        }

    def _measure_psf_moments(self, psf_obs):
        return self._measure_moments(
            psf_obs,
            self.psf_ampars,
            doround=False,
        )

    def _measure_obj_moments(self, obs):
        if False:
            import images
            images.multiview(obs[0].image)
            if 'q'==raw_input('hit a key: '):
                stop


        return self._measure_moments(
            obs,
            self.ampars,
        )

    def _measure_moments(self, obslist, ampars ,doround=True):
        """
        measure adaptive moments
        """
        #ampars=self['admom_pars']
        ntry=ampars.pop('ntry',4)

        use_ccen=ampars.get('use_canonical_center',False)
        if use_ccen:

            if isinstance(obslist,ngmix.ObsList):
                obs=obslist[0]
            else:
                obs=obslist
     
            ccen=(numpy.array(obs.image.shape)-1.0)/2.0
            jold=obs.jacobian
            obs.jacobian = ngmix.Jacobian(
                row=ccen[0],
                col=ccen[1],
                dvdrow=jold.dvdrow,
                dudrow=jold.dudrow,
                dvdcol=jold.dvdcol,
                dudcol=jold.dudcol,

            )


        fitter=ngmix.admom.Admom(
            obslist,
            rng=self.rng,
            **ampars
        )

        for i in xrange(ntry):
            Tguess=self._get_guess(obslist)

            fitter.go(Tguess)

            res=fitter.get_result()

            if doround:
                if res['flags'] != 0:
                    continue
                self._set_am_flux(obslist,fitter)

            if res['flags'] == 0:
                break

        if use_ccen:
            obs.jacobian=jold


        if res['flags'] != 0:
            raise TryAgainError("        admom failed")

        self._set_some_pars(res)
        #if doround:
        #    self._set_round_s2n(obslist,fitter)
        return res, fitter

    def _get_guess(self, obslist):
        if isinstance(obslist,ngmix.observation.ObsList):
            j=obslist[0].jacobian
        else:
            j=obslist.jacobian

        scale=j.get_scale()
        Tguess = 10.0*scale**2
        return Tguess
        
    def _get_guess_old(self, obs):

        rng=self.rng

        scale=obs.jacobian.get_scale()
        Tguess = 10.0*scale**2
        #Tguess = 4.0*scale**2

        pars=zeros(6)
        pars[0:0+2] = rng.uniform(low=-0.5*scale, high=0.5*scale, size=2)
        pars[2:2+2] = rng.uniform(low=-0.3, high=0.3, size=2)
        pars[4]     = Tguess*(1.0 + rng.uniform(low=-0.1, high=0.1))
        pars[5]     = 1.0

        guess=ngmix.GMixModel(pars, "gauss")
        return guess


    def _get_metacal(self, obslist):
        mcpars=self['metacal_pars']
        odict=ngmix.metacal.get_all_metacal(
            obslist,
            rng=self.rng,
            **mcpars
        )

        return odict

    def _do_metacal(self, odict):
        res={}

        for type in self.metacal_types:
            obslist=odict[type]

            tres,fitter=self._measure_obj_moments(obslist)
            if tres['flags'] != 0:
                raise TryAgainError("        bad T")

            if 'fit_model' in self:
                self._set_template_flux(obslist, tres)

            if type=='noshear':
                pres = self._measure_psfs(obslist)
                tres['psfrec_g'] = pres['g']
                tres['psfrec_T'] = pres['T']

                if False:
                    import images
                    images.multiview(obslist[0].image)
                    if 'q'==raw_input('hit a key: '):
                        stop

            res[type]=tres

        return res

    def _fit_flux_model(self, obs):

        mconf=self['max_pars']

        runner=self._get_runner(obs)

        try:

            runner.go(ntry=mconf['ntry'])

            fitter=runner.get_fitter() 
            res=fitter.get_result()
            if res['flags'] != 0:
                raise TryAgainError("failed to fit galaxy")

        except ngmix.GMixRangeError as err:
            raise TryAgainError("failed to fit galaxy: %s" % str(err))

        return fitter


    def _get_runner(self, obs):

        guesser=ngmix.guessers.PriorGuesser(self.prior)

        mconf=self['max_pars']
        runner=ngmix.galsimfit.GalsimRunner(
            obs,
            self['fit_model'],
            guesser,
            lm_pars=mconf['lm_pars'],
            prior=self.prior,
        )
        return runner

    def _set_template_flux(self, obs, res):

        try:
            # we use simulate_err=True since there is correlated
            # noise in the metacal images

            # How would we do this for multi-epoch?

            j=obs.jacobian
            scale = j.get_scale()
            jrow, jcol = j.get_cen()
            drow, dcol = res['pars'][0:2]

            row, col = jrow+drow/scale, jcol+dcol/scale

            # for the galsim object, we need to understand what this shift
            # means relative to the canonical galsim center

            # galsim true center in pixels
            gcen = (numpy.array(obs.image.shape)-1.0)/2.0

            offset_pixels = row-gcen[0], col-gcen[1]

            drowsky = offset_pixels[0]*scale
            dcolsky = offset_pixels[1]*scale

            model = self.model_galsim_obj.shift(dx=dcolsky, dy=drowsky)
            ffitter=ngmix.galsimfit.GalsimTemplateFluxFitter(
                obs,
                model,
                obs.psf.galsim_obj,
                simulate_err=True,
                rng=self.rng,
            )
            ffitter.go()
            fres=ffitter.get_result()
            if fres['flags'] != 0:
                raise TryAgainError("        template flux failure")

            res['flux'] = fres['flux']
            res['flux_err'] = fres['flux_err']
            res['flux_s2n'] = fres['flux']/fres['flux_err']
        except ngmix.GMixRangeError as err:
            raise TryAgainError(str(err))

 
    def _set_am_flux(self, obslist, amfitter):
        try:
            gmix=amfitter.get_gmix()

            for obs in obslist:
                obs.set_gmix(gmix)

            # we use simulate_err=True since there is correlated
            # noise in the metacal images
            fitter=ngmix.fitting.TemplateFluxFitter(
                obslist,
                rng=self.rng,
                simulate_err=True,
            )
            fitter.go()

            fres=fitter.get_result()
            if fres['flags'] != 0:
                res['flags'] = fres
                raise TryAgainError("could not get flux")

            res=amfitter.get_result()
            res['am_flux']=fres['flux']
            res['am_flux_err']=fres['flux_err']
            res['am_flux_s2n']=fres['flux']/fres['flux_err']

        except ngmix.GMixRangeError as err:
            raise TryAgainError(str(err))

    '''
    def _set_round_s2n(self, obslist, fitter):

        try:
            gm  = fitter.get_gmix()
            gmr = gm.make_round()

            e1r,e2r,T_r=gmr.get_e1e2T()

            res=fitter.get_result()
            flux=res['am_flux']
            gmr.set_flux(flux)

            s2n_sum=0.0
            for obs in obslist:
                s2n_sum += gmr.get_model_s2n(obs)

            if s2n_sum < 0.0:
                s2n_sum = 0.0

            res['s2n_r']=numpy.sqrt(s2n_sum)
            res['T_r'] = T_r

        except ngmix.GMixRangeError as err:
            raise TryAgainError(str(err))
    '''

    def _set_some_pars(self, res):
        res['g']     = res['e']
        res['g_cov'] = res['e_cov']

    def _copy_to_output(self, res, i):
        """
        copy parameters specific to this class
        """

        super(SimpleFitterBase,self)._copy_to_output(res, i)

        d=self.data

        pres=res['psf']
        d['psfrec_g'][i] = pres['g']
        d['psfrec_T'][i] = pres['T']

        if 'fit_model' in self:
            n=self._get_namer()
            model_res=res['model_result']
            for name in ['pars','sums_cov','g','g_cov','s2n_r']:
                d[n(name)][i] = model_res[name]


        for type in self.metacal_types:

            tres=res[type]
            if type=='noshear':
                back=''
            else:
                back='_%s' % type

            d['mcal_pars%s' % back][i] = tres['pars']
            d['mcal_g%s' % back][i] = tres['g']
            d['mcal_g_cov%s' % back][i] = tres['g_cov']

            if 'flux' in tres:
                d['mcal_flux%s' % back][i] = tres['flux']
                d['mcal_flux_s2n%s' % back][i] = tres['flux_s2n']

            d['mcal_am_flux%s' % back][i] = tres['am_flux']
            d['mcal_am_flux_s2n%s' % back][i] = tres['am_flux_s2n']


            d['mcal_am_flux%s' % back][i] = tres['am_flux']
            d['mcal_am_flux_s2n%s' % back][i] = tres['am_flux_s2n']

            d['mcal_s2n%s' % back][i] = tres['s2n']

            if 'T_r' in tres:
                d['mcal_s2n_r%s' % back][i] = tres['s2n_r']
                d['mcal_T_r%s' % back][i] = tres['T_r']

            d['mcal_numiter%s' % back][i] = tres['numiter']

            if type=='noshear':
                for p in ['sums_cov','wsum','psfrec_g','psfrec_T']:

                    if p in tres:
                        name='mcal_%s' % p
                        d[name][i] = tres[p]


    def _get_dtype(self):
        """
        get the dtype for the output struct
        """

        # always 6 for am
        npars=6

        dt=super(MetacalMomentsAM,self)._get_dtype()
        dt += [
            ('psfrec_g','f8',2),
            ('psfrec_T','f8'),
        ]

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
                    ('mcal_wsum','f8'),
                    ('mcal_sums_cov','f8',(npars,npars)),
                    ('mcal_psfrec_g','f8',2),
                    ('mcal_psfrec_T','f8'),
                ]

            dt += [
                ('mcal_s2n%s' % back,'f8'),
                #('mcal_s2n_r%s' % back,'f8'),
                #('mcal_T_r%s' % back,'f8'),

                ('mcal_am_flux%s' % back,'f8'),
                ('mcal_am_flux_s2n%s' % back,'f8'),
            ]

            if 'fit_model' in self:
                dt += [
                    ('mcal_flux%s' % back,'f8'),
                    ('mcal_flux_s2n%s' % back,'f8'),
                ]

            dt += [('mcal_numiter%s' % back,'i4')]

        return dt

    def _print_res(self,res):
        """
        print some stats
        """

        subres=res['noshear']

        logger.debug("")

        if 'prefit' in res:
            preres=res['prefit']
            logger.debug("    s2n: %g" % preres['s2n'])

        if 'fit_model' in self:
            mres=res['model_result']
            pars=mres['pars']
            perr=mres['pars_err']
            log_pars(pars,front='    model:')
            log_pars(perr,front='          ')

        s2n=subres['s2n']
        g=subres['g']

        cov = subres['g_cov'].clip(min=0, max=None)
        gerr=diag(sqrt(cov))

        logger.debug('    true r50: %(r50_true)g flux: %(flux_true)g s2n: %(s2n_true)g' % res)
        logger.debug("    mcal s2n: %g  e:  %g +/- %g  %g +/- %g" % (s2n,g[0],gerr[0],g[1],gerr[1]))
        if 'flux' in subres:
            logger.debug("        flux: %g +/- %g flux_s2n:  %g" % (subres['flux'],subres['flux_err'],subres['flux_s2n']))
        logger.debug("     am numiter: %d flux: %g +/- %g flux_s2n:  %g" % (subres['numiter'], subres['am_flux'],subres['am_flux_err'],subres['am_flux_s2n']))


    def _do_plots(self, obs, gmix=None):
        import images

        if gmix is not None:
            gm=gmix.copy()
            gm.set_psum(obs.image.sum())

            model=gm.make_image(obs.image.shape, jacobian=obs.jacobian)
            images.compare_images(
                obs.image,
                model,
                label1='image',
                label2='model',
                file='/u/ki/esheldon/public_html/tmp/plots/tmp.png',
                dims=[1000,1000],
            )

        else:
            images.multiview(
                obs.image,
                title='image',
                file='/u/ki/esheldon/public_html/tmp/plots/tmp.png',
                dims=[1000,1000],
            )

        if 'q'==raw_input("hit a key: "):
            stop


class MetacalMomentsAMMOFSub(MetacalMomentsAM):
    def _dofit(self, allobs):
        import minimof

        try:
            mm=minimof.MiniMOF(
                self['mof'],
                allobs,
                rng=self.rng,
            )
            mm.go()
            mm_res = mm.get_result()
            if not mm_res['converged']:
                raise TryAgainError("MOF did not converge")
        except BootPSFFailure as err:
            raise TryAgainError("MOF psf failure: '%s'" % str(err))
        except BootGalFailure as err:
            raise TryAgainError("MOF gal failure: '%s'" % str(err))

        # assume first is the central
        corr_obslist = mm.get_corrected_obs(0)

        return super(MetacalMomentsAMMOFSub,self)._dofit(corr_obslist)


class AMFitter(MetacalMomentsAM):
    def _setup(self, *args, **kw):
        super(MetacalMomentsAM,self)._setup(*args, **kw)
        self._set_mompars()

        self['min_s2n'] = self.get('min_s2n',0.0)


    def _dofit(self, obslist):

        res,fitter=self._measure_moments(obslist, self.ampars)
        res['flags']=0
        return res

    def _dofit(self, obslist):

        psfres = self._measure_psfs(obslist)

        res,fitter=self._measure_moments(obslist, self.ampars)
        self._set_am_flux(obslist,fitter)
        res['flags']=0

        res['psf'] = psfres

        return res


    def _get_dtype(self):
        """
        get the dtype for the output struct
        """
        npars=self['npars']

        # super of super
        dt=super(SimpleFitterBase,self)._get_dtype()
        dt += [
            ('pars','f8',npars),
            ('sums_cov','f8',(npars,npars)),
            ('T','f8'),
            ('T_err','f8'),
            ('flux','f8'),
            ('flux_err','f8'),
            ('g','f8',2),
            ('g_cov','f8',(2,2)),
            ('s2n','f8'),
            ('numiter','i4'),
        ]

        return dt

    def _copy_to_output(self, res, i):
        """
        copy parameters specific to this class
        """

        # note copying super of our super, since
        # we didn't do a regular fit
        super(SimpleFitterBase,self)._copy_to_output(res, i)

        d=self.data

        ckeys=[
            'pars','sums_cov',
            'g','g_cov',
            's2n',
            'numiter',
        ]

        d['flux'][i] = res['am_flux']
        d['flux_err'][i] = res['am_flux_err']
        d['T'][i] = res['T']
        d['T_err'][i] = res['T_err']
        for key in ckeys:
            if key in res:
                d[key][i] = res[key]

    def _print_res(self,res):
        """
        print some stats
        """

        if 'flux_s2n' in res:
            logger.debug("    flux s2n: %g" % res['flux_s2n'])
        logger.debug("    e1e2:  %g %g" % tuple(res['g']))
        logger.debug("    e_err: %g" % numpy.sqrt(res['g_cov'][0,0]))

        log_pars(res['pars'],      front='        pars: ')

        logger.debug('        true r50: %(r50_true)g flux: %(flux_true)g s2n: %(s2n_true)g' % res)


def get_shears(step):
    return {
        'noshear':None,
        '1p':Shape( step, 0.0),
        '1m':Shape(-step, 0.0),
        '2p':Shape( 0.0,  step),
        '2m':Shape( 0.0, -step),
    }

class KMetacal(ngmix.metacal.Metacal):
    """
    currently only works with symmetrized psf

    """
    def __init__(self, obs, dim, dk, step=0.01, **kw):

        self.dk=dk
        self.dim=dim
        self.step=step
        self.shears=get_shears(step)
        self.dilation=1.0 + 2.0*step

        super(KMetacal,self).__init__(obs, **kw)

        assert self.symmetrize_psf

    def _set_data(self):
        """
        only with pixel for now
        """
        super(KMetacal,self)._set_data()

        self._dilated_psf_nopix = self.psf_int_nopix.dilate(self.dilation)
        self._dilated_psf = galsim.Convolve(self.psf_int_nopix,self.pixel)
        self._dilated_psf_kr, self._dilated_psf_ki = self._make_kimages(self._dilated_psf)

    def get_kobs_galshear(self, shtype):
        """
        This is the case where we shear the image, for calculating R

        parameters
        ----------
        shear: ngmix.Shape
            The shear to apply

        get_unsheared: bool
            Get an observation only convolved by the target psf, not
            sheared
        """

        kr, ki = self.get_target_kimages(shtype)
        newobs = self._make_kobs(kr, ki)

        return newobs


    def get_target_kimages(self, shtype):
        """
        get the target image, convolved with the specified psf
        and possibly sheared

        parameters
        ----------
        psf_obj: A galsim object
            psf object by which to convolve.  An interpolated image,
            or surface brightness profile
        shear: ngmix.Shape, optional
            The shear to apply

        returns
        -------
        galsim image object
        """

        shear=self.shears[shtype]
        imconv = self._get_target_gal_obj(shear=shear)

        # this should carry over the wcs
        kr,ki = self._make_kimages(imconv)
        return kr, ki

    def _get_target_gal_obj(self, shear=None):
        if shear is not None:
            shim_nopsf = self.get_sheared_image_nopsf(shear)
        else:
            shim_nopsf = self.image_int_nopsf

        imconv = galsim.Convolve([shim_nopsf, self._dilated_psf])

        return imconv

    def _make_kimages(self, obj):
        # this should carry over the wcs
        kr,ki = obj.drawKImage(
            dtype=numpy.float64,
            nx=self.dim,
            ny=self.dim,
            scale=self.dk,
        )

        return kr, ki

    def _make_kobs(self, kr, ki):
        from ngmix.observation import KObservation

        psf_kr, psf_ki = self._dilated_psf_kr.copy(), self._dilated_psf_ki.copy()

        weight = kr.copy()
        medweight = numpy.median(self.obs.weight)
        weight.array[:,:] = 0.5*medweight

        # parseval's theorem
        #weight *= (1.0/weight.array.size)
        weight *= (1.0/self.obs.image.size)


        psf_medweight = numpy.median(self.obs.psf.weight)
        psf_weight = psf_kr.copy()
        psf_weight.array[:,:] = 0.5*psf_medweight

        psf_kobs = KObservation(
            psf_kr,
            psf_ki,
            weight=psf_weight,
        )

        kobs = KObservation(
            kr,
            ki,
            weight=weight,
            psf=psf_kobs,
        )

        return kobs

class MetacalMomentsFixed(MetacalMomentsAM):
    """
    fixed weight function
    """

    def _set_mompars(self):
        super(MetacalMomentsFixed,self)._set_mompars()
        wpars=self['weight']

        self.weight=ngmix.GMixModel(
            [0.0, 0.0, 0.0, 0.0, wpars['T'], 1.0],
            wpars['model'],
        )

    def _measure_psf_moments(self, psf_obs):
        """
        want to run admom for psf
        """
        return super(MetacalMomentsFixed,self)._measure_moments(
            psf_obs,
            self.psf_ampars,
            doround=False,
        )


    def _measure_moments(self, obslist, junk, doround=True):
        """
        measure adaptive moments
        """

        if isinstance(obslist,ngmix.ObsList):
            obs=obslist[0]
        else:
            obs=obslist
            
        wpars=self['weight']
        if wpars['use_canonical_center']:
        
            ccen=(numpy.array(obs.image.shape)-1.0)/2.0
            jold=obs.jacobian
            obs.jacobian = ngmix.Jacobian(
                row=ccen[0],
                col=ccen[1],
                dvdrow=jold.dvdrow,
                dudrow=jold.dudrow,
                dvdcol=jold.dvdcol,
                dudcol=jold.dudcol,

            )

        res = self.weight.get_weighted_moments(
            obs,
        )

        if wpars['use_canonical_center']:
            obs.jacobian=jold

        if res['flags'] != 0:
            raise TryAgainError("        moments failed")

        res['e'] = res['pars'][2:2+2]
        res['e_cov'] = res['pars_cov'][2:2+2, 2:2+2]

        pars_cov=res['pars_cov']
        flux_sum=res['pars'][5]
        fvar_sum=pars_cov[5,5]

        if fvar_sum > 0.0:

            flux_err = numpy.sqrt(fvar_sum)
            res['s2n'] = flux_sum/flux_err

            # error on each shape component from BJ02 for gaussians
            # assumes round

            res['e_err_r'] = 2.0/res['s2n']
        else:
            res['flags'] = 0x40
            flux_err = 9999.0

        res['g']     = res['e']
        res['g_cov'] = res['e_cov']
        res['T'] = res['pars'][4]
        res['am_flux'] = flux_sum
        res['am_flux_err'] = flux_err
        res['am_flux_s2n'] = res['s2n']
        res['numiter']=1

        fitter=None
        return res, fitter


