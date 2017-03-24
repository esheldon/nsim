"""
moment based code
"""
from __future__ import print_function
try:
    xrange
except:
    xrange=range

import numpy
from numpy import array, zeros, ones, log, log10, exp, sqrt, diag
from numpy import where, isfinite

import ngmix
from ngmix.fitting import print_pars
from ngmix.shape import Shape

from . import util
from .util import TryAgainError

from .fitters import SimpleFitterBase

import galsim


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
        print("doing types:",self.metacal_types)


        self['min_s2n'] = self.get('min_s2n',0.0)


    def _dofit(self, imdict):

        obs=imdict['obs']

        psfres, fitter = self._measure_moments(obs.psf, doround=False)

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

        obsdict=self._get_metacal(obs)

        res=self._do_metacal(obsdict)

        res['psf'] = {}
        res['psf']['g'] = psfres['g']
        res['psf']['T'] = psfres['T']

        if 'fit_model' in self:
            res['model_result'] = flux_result

        res['flags']=0
        return res

    def _measure_moments(self, obs, doround=True):
        """
        measure adaptive moments
        """
        ampars=self['admom_pars']
        ntry=ampars.pop('ntry',4)

        fitter=ngmix.admom.Admom(obs, **ampars)

        for i in xrange(ntry):
            guess=self._get_guess(obs)

            fitter.go(guess)

            res=fitter.get_result()

            if doround:
                if res['flags'] != 0:
                    continue
                self._set_am_flux(obs,fitter)

            if res['flags'] == 0:
                break

        if res['flags'] != 0:
            raise TryAgainError("        admom failed")


        self._set_some_pars(res)
        if doround:
            self._set_round_s2n(obs,fitter)
        return res, fitter

    def _get_guess(self, obs):

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


    def _get_metacal(self, obs):
        mcpars=self['metacal_pars']
        odict=ngmix.metacal.get_all_metacal(
            obs,
            rng=self.rng,
            **mcpars
        )

        return odict

    def _do_metacal(self, odict):
        res={}

        for type in self.metacal_types:
            #print("    doing metacal:",type)
            obs=odict[type]

            tres,fitter=self._measure_moments(obs)
            if tres['flags'] != 0:
                raise TryAgainError("        bad T")

            self._set_template_flux(obs, tres)

            if type=='noshear':
                pres,fitter=self._measure_moments(obs.psf, doround=False)
                tres['psfrec_g'] = pres['g']
                tres['psfrec_T'] = pres['T']

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

 
    def _set_am_flux(self, obs, amfitter):
        try:
            gmix=amfitter.get_gmix()

            obs.set_gmix(gmix)

            # we use simulate_err=True since there is correlated
            # noise in the metacal images
            fitter=ngmix.fitting.TemplateFluxFitter(obs, simulate_err=True)
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

    def _set_round_s2n(self, obs, fitter):

        try:
            gm  = fitter.get_gmix()
            gmr = gm.make_round()

            e1r,e2r,T_r=gmr.get_e1e2T()

            res=fitter.get_result()
            flux=res['am_flux']
            gmr.set_flux(flux)

            res['s2n_r']=gmr.get_model_s2n(obs)
            res['T_r'] = T_r

        except ngmix.GMixRangeError as err:
            raise TryAgainError(str(err))

    def _set_some_pars(self, res):
        res['g']     = res['e']
        res['g_cov'] = res['e_cov']

        # not right pars cov
        res['pars_cov']=res['sums_cov']*0 + 9999.e9

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
            for name in ['pars','pars_cov','g','g_cov','s2n_r']:
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

            d['mcal_flux%s' % back][i] = tres['flux']
            d['mcal_flux_s2n%s' % back][i] = tres['flux_s2n']

            d['mcal_am_flux%s' % back][i] = tres['am_flux']
            d['mcal_am_flux_s2n%s' % back][i] = tres['am_flux_s2n']

            d['mcal_s2n%s' % back][i] = tres['s2n']
            d['mcal_s2n_r%s' % back][i] = tres['s2n_r']
            d['mcal_T_r%s' % back][i] = tres['T_r']

            d['mcal_numiter%s' % back][i] = tres['numiter']

            if type=='noshear':
                for p in ['pars_cov','wsum','psfrec_g','psfrec_T']:

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
                    ('mcal_pars_cov','f8',(npars,npars)),
                    ('mcal_psfrec_g','f8',2),
                    ('mcal_psfrec_T','f8'),
                ]

            dt += [
                ('mcal_s2n%s' % back,'f8'),
                ('mcal_s2n_r%s' % back,'f8'),
                ('mcal_T_r%s' % back,'f8'),
                ('mcal_flux%s' % back,'f8'),
                ('mcal_flux_s2n%s' % back,'f8'),

                ('mcal_am_flux%s' % back,'f8'),
                ('mcal_am_flux_s2n%s' % back,'f8'),
            ]

            dt += [('mcal_numiter%s' % back,'i4')]

        return dt

    def _print_res(self,res):
        """
        print some stats
        """

        subres=res['noshear']

        print()

        if 'prefit' in res:
            preres=res['prefit']
            print("    s2n: %g" % preres['s2n'])

        if 'fit_model' in self:
            mres=res['model_result']
            pars=mres['pars']
            perr=mres['pars_err']
            print_pars(pars,front='    model:')
            print_pars(perr,front='          ')

        s2n=subres['s2n']
        g=subres['g']

        cov = subres['g_cov'].clip(min=0, max=None)
        gerr=diag(sqrt(cov))

        print("    mcal s2n: %g  e:  %g +/- %g  %g +/- %g" % (s2n,g[0],gerr[0],g[1],gerr[1]))
        print("        flux: %g +/- %g flux_s2n:  %g" % (subres['flux'],subres['flux_err'],subres['flux_s2n']))
        print("     am flux: %g +/- %g flux_s2n:  %g" % (subres['am_flux'],subres['am_flux_err'],subres['am_flux_s2n']))


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


class AMFitter(MetacalMomentsAM):
    def _dofit(self, imdict):

        obs=imdict['obs']

        res,fitter=self._measure_moments(obs)
        res['flags']=0
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
            ('pars_cov','f8',(npars,npars)),
            ('flux','f8'),
            ('flux_s2n','f8'),
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
            'pars','pars_cov',
            'flux','flux_s2n',
            'g','g_cov',
            's2n',
            'numiter',
        ]

        for key in ckeys:
            d[key][i] = res[key]

    def _print_res(self,res):
        """
        print some stats
        """

        print("    flux s2n: %g" % res['flux_s2n'])
        print("    e1e2:  %g %g" % tuple(res['g']))
        print("    e_err: %g" % numpy.sqrt(res['g_cov'][0,0]))

        print_pars(res['pars'],      front='        pars: ')

        print('        true: ', res['pars_true'])


#
# failed experiments
#

class MetacalMomentsFixed(SimpleFitterBase):
    """
    - might be effects not captured by the galsim
      Convolve
    - probably should dilate while we are symmetrizing
    """
    def _setup(self, *args, **kw):
        super(MetacalMomentsFixed,self)._setup(*args, **kw)

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
        print("doing types:",self.metacal_types)


        self['min_s2n'] = self.get('min_s2n',0.0)

    def _dofit(self, imdict):

        obs=imdict['obs']

        self._set_weight(obs) 

        odict=self._get_metacal(obs)

        #self._set_F(self.kmcal.dim)

        res=self._do_metacal(odict)
        res['flags']=0

        return res

    def _do_metacal(self, odict):
        res={}

        for type in self.metacal_types:
            #print("    doing metacal:",type)
            obs=odict[type]

            tres,fitter=self._measure_moments(obs)
            if tres['flags'] != 0:
                raise TryAgainError("        bad T")


            if type=='noshear':
                #pres,fitter=self._measure_moments(obs.psf_nopix)
                pres,fitter=self._measure_moments(obs.psf, doround=False)
                tres['psfrec_g'] = pres['g']
                tres['psfrec_T'] = pres['T']

            res[type]=tres

        return res

    def _measure_moments(self, obs):
        """
        moments of the power spectrum

        todo: subtract noise ps. Deal with errors
        """

        jac=obs.jacobian
        cen=jac.get_cen()
        dims=obs.image.shape

        scale=jac.get_scale()

        min_rdiff=min(cen[0], (dims[0]-1)-cen[0] )
        min_cdiff=min(cen[1], (dims[1]-1)-cen[1] )
        min_diff=min(min_rdiff, min_cdiff)

        rmax=min_diff*scale

        # flags would only be set if there were ivar <= 0
        res=self.wt_gmix.get_weighted_moments(obs, rmax=rmax)

        if False:
            momslow=self._get_moments_slow(obs)
            print_pars(res['pars'],  front="    pars fast:")
            print_pars(momslow,      front="    pars slow:")
            print_pars(momslow-res['pars'], front="    diff:     ",fmt='%.16g')

        if False:
            self._do_plots(obs)

        res['pars_err']=sqrt(diag(res['pars_cov']))

        e=zeros(2) - 9999.0
        ecov=zeros( (2,2)) + 9999.0
        flux=-9999.0
        flux_s2n=-9999.0
        T=-9999.0

        pars=res['pars']
        pars_cov=res['pars_cov']
        if pars[4] != 0 and pars[5] != 0.0:
            e[0]=pars[2]/pars[4]
            e[1]=pars[3]/pars[4]
            T = pars[4]/pars[5]
        else:
            res['flags']=2**0

        if pars_cov[5,5] != 0.0:
            flux     = pars[5]/res['wsum']
            flux_s2n = pars[5]/sqrt(pars_cov[5,5])
        else:
            res['flags']=2**1
            
        res['g']        = e
        res['g_cov']    = ecov
        res['flux']     = flux
        res['flux_s2n'] = flux_s2n
        res['T']        = T

        return res, None

    def _set_weight(self, obs):
        """
        for now used fixed gaussian weight function
        """
        self.wt_gmix = _get_weight(obs, self)

    def _get_dtype(self):
        """
        get the dtype for the output struct
        """
        npars=self['npars']

        # super of super
        dt=super(MetacalMomentsFixed,self)._get_dtype()

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
                    ('mcal_pars_cov','f8',(npars,npars)),
                    ('mcal_psfrec_g','f8',2),
                    ('mcal_psfrec_T','f8'),
                ]

            dt += [
                ('mcal_s2n%s' % back,'f8'),
                ('mcal_s2n_r%s' % back,'f8'),
                ('mcal_T_r%s' % back,'f8'),
                ('mcal_flux%s' % back,'f8'),
                ('mcal_flux_s2n%s' % back,'f8'),
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

        for type in self.metacal_types:

            tres=res[type]
            if type=='noshear':
                back=''
            else:
                back='_%s' % type

            d['mcal_pars%s' % back][i] = tres['pars']
            d['mcal_g%s' % back][i] = tres['g']
            d['mcal_g_cov%s' % back][i] = tres['g_cov']

            d['mcal_flux%s' % back][i] = tres['flux']
            d['mcal_flux_s2n%s' % back][i] = tres['flux_s2n']

            d['mcal_s2n%s' % back][i] = tres['s2n']
            d['mcal_s2n_r%s' % back][i] = tres['s2n_r']
            d['mcal_T_r%s' % back][i] = tres['T_r']

            if type=='noshear':
                for p in ['pars_cov','wsum','psfrec_g','psfrec_T']:

                    if p in tres:
                        name='mcal_%s' % p
                        d[name][i] = tres[p]

    def _print_res(self,res):
        """
        print some stats
        """

        subres=res['noshear']

        print()

        if 'prefit' in res:
            preres=res['prefit']
            print("    s2n: %g" % preres['s2n'])

        s2n=subres['s2n']
        g=subres['g']

        cov = subres['g_cov'].clip(min=0, max=None)
        gerr=diag(sqrt(cov))

        print("    mcal s2n: %g  e:  %g +/- %g  %g +/- %g" % (s2n,g[0],gerr[0],g[1],gerr[1]))
        #print_pars(subres['pars'],      front='        pars: ')



    def _do_plots(self, obs):
        import biggles
        import images

        tab=biggles.Table(2, 1)


        wtim=self.wt_gmix.make_image(obs.image.shape, jacobian=obs.jacobian)

        tab[0,0]=images.multiview(
            obs.image,
            title='image',
            show=False,
        )
        tab[1,0]=images.multiview(
            obs.image*wtim,
            title='weighted image',
            show=False,
        )

        tab.write_img(1000,1500,'/u/ki/esheldon/public_html/tmp/plots/tmp.png')
        if 'q'==raw_input("hit a key: "):
            stop

    def _get_moments_slow(self, obs):
        cen=obs.jacobian.get_cen()
        dims=obs.image.shape
        rows,cols=numpy.mgrid[
            0:dims[0],
            0:dims[1],
        ]

        rows=rows.astype('f8')-cen[0]
        cols=cols.astype('f8')-cen[1]

        self.F0 = rows.copy()
        self.F1 = cols.copy()
        self.F2 = cols**2 - rows**2
        self.F3 = 2.0*cols*rows
        self.F4 = cols**2 + rows**2
        self.F5 = rows*0 + 1

        wt=self.wt_gmix.make_image(obs.image.shape, jacobian=obs.jacobian)

        wtim=wt*obs.image

        pars=zeros(6)
        pars[0] = (self.F0*wtim).sum()
        pars[1] = (self.F1*wtim).sum()
        pars[2] = (self.F2*wtim).sum()
        pars[3] = (self.F3*wtim).sum()
        pars[4] = (self.F4*wtim).sum()
        pars[5] = wtim.sum()

 
        return pars

    def _get_guess(self, obs):


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


class MetacalGaussK(MetacalMomentsAM):
    def _dofit(self, imdict):

        obs=imdict['obs']
        if False:
            self._do_plots(obs)

        self.flux_guess = obs.image.sum()
        self.psf_flux_guess = obs.psf.image.sum()

        obsdict=self._get_metacal(obs)

        res=self._do_metacal(obsdict)

        res['flags']=0
        return res

    def _get_metacal(self, obs):
        """
        simplified version of fixnoise, using minus shear
        """
        noise_obs = ngmix.simobs.simulate_obs(None, obs)

        # just doing this to get the dim, dk.  This is wasting time
        # but I'm not sure how much
        mb_iilist, dim, dk = ngmix.observation.make_iilist(obs)

        odict={}

        mcal_pars=self['metacal_pars']
        km = KMetacal(obs,  dim, dk, **mcal_pars)
        kmnoise = KMetacal(noise_obs,  dim, dk, **mcal_pars)

        for type in self.metacal_types:

            kobs = km.get_kobs_galshear(type)

            if type=='1p':
                nkobs = kmnoise.get_kobs_galshear("1m")
            elif type=="1m":
                nkobs = kmnoise.get_kobs_galshear("1p")
            elif type=='2p':
                nkobs = kmnoise.get_kobs_galshear("2m")
            elif type=="2m":
                nkobs = kmnoise.get_kobs_galshear("2p")
            elif type=="noshear":
                nkobs = kmnoise.get_kobs_galshear("noshear")
            else:
                raise ValueError("bad type: '%s'" % type)


            kobs.kr.array[:,:] += nkobs.kr.array[:,:]
            kobs.ki.array[:,:] += nkobs.ki.array[:,:]

            # doubling the variance
            kobs.weight *= 0.5

            #var = 0.5*(nkobs.kr.array.var() + nkobs.ki.array.var())
            #kobs.weight.array[:,:] = 1.0/var

            odict[type] = kobs

        return odict

    def _do_metacal(self, odict):
        res={}

        for type in self.metacal_types:
            obs=odict[type]

            tres,fitter=self._fit_one(obs, self.flux_guess)
            if tres['flags'] != 0:
                raise TryAgainError("        failed to fit")

            if type=='noshear':
                pres,fitter=self._fit_one(obs.psf, self.psf_flux_guess)
                tres['psfrec_g'] = pres['g']
                tres['psfrec_T'] = pres['pars'][4]

            res[type]=tres

        return res


    def _fit_one(self, obs, flux_guess):
        """
        actually do a max like fit, but we want to share
        this simpler interface
        """
        mpars=self['max_pars']
        lmpars=mpars['lmpars']

        ntry=mpars.pop('ntry',4)

        fitter=ngmix.fitting.LMGaussK(obs, lmpars=lmpars)

        for i in xrange(ntry):
            guess=self._get_guess(obs, flux_guess)

            fitter.go(guess)

            res=fitter.get_result()

            if res['flags'] == 0:
                break

        if res['flags'] != 0:
            raise TryAgainError("        lm failed")

        pcov=res['pars_cov']
        if pcov[5,5] > 0:
            pars=res['pars']
            res['flux'] = pars[5]
            res['flux_s2n'] = pars[5]/pcov[5,5]
        else:
            raise TryAgainError("        bad flux error")

        res['s2n'] = res['s2n_w']
        res['s2n_r'] = res['s2n']
        res['numiter'] = res['nfev']
        res['T_r'] = -9999.0

        return res, fitter

    def _get_guess(self, obs, flux_guess):


        rng=self.rng

        scale=obs.jacobian.get_scale()
        Tguess = 10.0*scale**2
        #Tguess = 4.0*scale**2

        pars=zeros(6)
        pars[0:0+2] = rng.uniform(low=-0.5*scale, high=0.5*scale, size=2)
        pars[2:2+2] = rng.uniform(low=-0.3, high=0.3, size=2)
        pars[4]     = Tguess*(1.0 + rng.uniform(low=-0.1, high=0.1))
        pars[5]     = flux_guess*(1.0 + rng.uniform(low=-0.1, high=0.1))

        return pars



class MetacalMomentsAMFixed(MetacalMomentsAM):
    def _do_metacal(self, odict):
        res={}

        self._set_weight(self.pre_fitter)

        for type,obs in odict.iteritems():

            tres=self._measure_moments_fix(obs)
            if tres['flags'] != 0:
                raise TryAgainError("        bad T")


            if type=='noshear':
                pres,fitter=self._measure_moments(obs.psf)
                tres['psfrec_g'] = pres['g']
                tres['psfrec_T'] = pres['T']

            res[type]=tres

        return res

    def _measure_moments_fix(self, obs):
        """
        measure deweighted moments
        """

        try:
            res=self.wt_gmix.get_weighted_moments(obs)
        except ngmix.GMixRangeError:
            raise TryAgainError("bad weight")

        # these are not normalized
        wpars=res['pars']
        M1 = wpars[2]/wpars[5]
        M2 = wpars[3]/wpars[5]
        T  = wpars[4]/wpars[5]

        Minv = self._get_Minv(M1, M2, T)

        Minv_deweight = Minv - self.Mwt_inv

        try:
            M = numpy.linalg.inv(Minv_deweight )
        except numpy.linalg.LinAlgError:
            raise TryAgainError("        could not invert deweighted moment matrix")

        Irr = M[0,0]
        Irc = M[0,1]
        Icc = M[1,1]

        M1 = Icc - Irr
        M2 = 2*Irc
        T = Icc + Irr

        if T <= 0.0:
            raise TryAgainError("        admom deweight gave T <= 0")

        # pack into similar structure from admom
        res['sums'] = res['pars']
        res['sums_cov'] = res['pars_cov']

        # these would normally represent the weight
        res['pars'] = res['sums'] * 0 - 9999.0

        cen=self.wt_gmix.get_cen()
        res['pars'][0:0+2] = cen
        res['pars'][2] = M1/T
        res['pars'][3] = M2/T
        res['pars'][4] = T
        res['pars'][5] = 1.0
        res['numiter'] = 1

        res=ngmix.admom.copy_result(res)
        self._set_some_pars(res)

        return res

    def _set_weight(self, fitter):
        """
        set the weight gaussian mixture, as well as 
        its inverse covariance matrix
        """

        self.wt_gmix = fitter.get_gmix()

        res=fitter.get_result()
        pars=res['pars']
        M1 = pars[2]
        M2 = pars[3]
        T  = pars[4]

        self.Mwt_inv = self._get_Minv(M1, M2, T)

    def _get_Minv(self, M1, M2, T):

        M = _M1M2T_to_matrix(M1, M2, T)
        try:
            Minv = numpy.linalg.pinv(M)
        except numpy.linalg.LinAlgError:
            raise TryAgainError("        could not invert moment matrix")

        return Minv


class MetacalMomentsAMMulti(MetacalMomentsAM):

    def _dofit(self, imdict):

        obs=imdict['obs']

        odict_list=self._get_metacal(obs)

        avgkeys=[
            'g','g_cov',
            'flux','flux_s2n',
            's2n','pars','pars_cov',
        ]
        for i,odict in enumerate(odict_list):
            tres=self._do_metacal(odict)

            if i==0:
                res=tres
            else:
                for type in res:
                    tsubres=tres[type]
                    subres=res[type]
                    
                    for key in avgkeys:
                        subres[key] += tsubres[key]

        nrand=len(odict_list)

        cov_fac=(1.0 + 1.0/nrand)/2.0
        s2n_fac = sqrt(1.0/cov_fac)

        for type in res:
            subres=res[type]
            for key in avgkeys:
                subres[key] /= nrand

            subres['s2n'] *= s2n_fac
            subres['flux_s2n'] *= s2n_fac

            subres['g_cov'] *= cov_fac
            subres['pars_cov'] *= cov_fac

        res['flags']=0

        return res

 
    def _get_metacal(self, obs):
        mcpars=self['metacal_pars']
        odict_list=get_all_metacal_multi(
            obs,
            rng=self.rng,
            **mcpars
        )

        return odict_list


class MetacalMomentsDeweight(MetacalMomentsFixed):

    def _measure_psf_am(self, obs):
        """
        moments of the power spectrum

        todo: subtract noise ps. Deal with errors
        """

        ampars=self['admom_pars']
        ntry=ampars['ntry']

        fitter=ngmix.admom.Admom(obs.psf, **ampars)

        for i in xrange(ntry):
            guess=self._get_guess(obs.psf)

            fitter.go(guess)
            res=fitter.get_result()

            if res['flags'] == 0:
                break

        if res['flags'] != 0:
            raise TryAgainError("        admom failed for psf")

        return res

    def _measure_moments(self, obs):
        """
        moments of the power spectrum

        todo: subtract noise ps. Deal with errors
        """

        # [ [Irr, Irc]
        #   [Irc, Icc] ]
        psf_res = self._measure_psf_am(obs)

        ppars=psf_res['pars']
        Mpsf = _M1M2T_to_matrix(ppars[2], ppars[3], ppars[4])

        res=super(MetacalMomentsDeweight,self)._measure_moments(obs)
        if res['flags'] != 0:
            raise TryAgainError("        bad moments")

        wpars=res['pars']
        M1 = wpars[2]/wpars[5]
        M2 = wpars[3]/wpars[5]
        T  = wpars[4]/wpars[5]

        Mmeas = _M1M2T_to_matrix(M1, M2, T)

        try:
            Mmeas_inv = numpy.linalg.pinv(Mmeas)
        except numpy.linalg.LinAlgError:
            raise TryAgainError("        could not invert observed moment matrix")


        Minv_deweight = Mmeas_inv - self.Mwt_inv

        try:
            Mdeweight = numpy.linalg.pinv(Minv_deweight )
        except numpy.linalg.LinAlgError:
            raise TryAgainError("        could not invert deweighted moment matrix")
        
        M = Mdeweight - Mpsf

        Irr = M[0,0]
        Irc = M[0,1]
        Icc = M[1,1]

        res['T'] = Irr + Icc
        if res['T']==0.0:
            res['flags'] = 2**0
        else:
            res['g'][0] = (Icc - Irr)/res['T']
            res['g'][1] = 2.0*Irc/res['T']

        return res, None


    def _set_weight(self, obs):
        """
        for now used fixed gaussian weight function
        """
        super(MetacalMomentsDeweight,self)._set_weight(obs)

        g1,g2,T = self.wt_gmix.get_g1g2T()

        M=zeros( (2,2) )
        Minv=zeros( (2,2) )
        M[0,0] = T/2.0
        M[1,1] = T/2.0

        Minv[0,0] = 1.0/M[0,0]
        Minv[1,1] = 1.0/M[1,1]

        self.Mwt = M

        self.Mwt_inv = Minv

    def _copy_to_output(self, res, i):
        """
        copy parameters specific to this class
        """

        # note copying super of our super, since
        # we didn't do a regular fit
        super(MetacalMomentsDeweight,self)._copy_to_output(res, i)

        d=self.data

        for type in self.metacal_types:
            if type=='noshear':
                back=''
            else:
                back='_%s' % type

            tres=res[type]
            d['mcal_T%s' % back][i] = tres['T']



    def _get_dtype(self):
        """
        get the dtype for the output struct
        """
        # super of super
        dt=super(MetacalMomentsDeweight,self)._get_dtype()

        for type in self.metacal_types:
            if type=='noshear':
                back=''
            else:
                back='_%s' % type

            dt += [('mcal_T%s' % back,'f8')]

        return dt

class NullerBase(object):
    def __init__(self, obs, **kw):
        self.obs=obs
        self.jacobian=self.obs.jacobian.copy()

        self.kw=kw

        self.npars=4


    def get_result(self):
        """
        get the result dictionary
        """
        if not hasattr(self,'_result'):
            raise RuntimeError("run go() first")
        
        return self._result

    def go(self, guess):
        """
        run the root finder and set the result
        """
        import scipy.optimize

        assert len(guess)==self.npars

        root_res=scipy.optimize.root(
            self.get_moments,
            guess,
            options=self.kw,
            method='lm',
        )

        if not root_res.success:
            self._result={'flags':1}
            return

        fit_pars=root_res.x

        res=self.get_moments(fit_pars, more=True)
        if res['s2n_denom_sum'] <= 0:
            self._result={'flags':1}
            return


        # rename some things
        res['sums'] = res['pars']
        res['sums_cov'] = res['pars_cov']
        res['pars'] = fit_pars

        if self.npars==2:
            res['g'] = -fit_pars[0:0+2]
        else:
            res['g'] = -fit_pars[2:2+2]

        # now some noise calculations

        res['s2n'] = res['s2n_numer_sum']/sqrt(res['s2n_denom_sum'])
        res['e_err'] = 2.0/res['s2n']

        gerr=0.5*res['e_err']
        res['g_cov'] = diag([ gerr**2 ]*2)

        res['flags'] = 0
        res['nfev'] = root_res.nfev

        self._result = res

class NullerGauss(NullerBase):
    def __init__(self, obs, T, **kw):

        super(NullerGauss,self).__init__(obs, **kw)
        self.T = T

        wtpars = [
            0.0,
            0.0,
            0.0,
            0.0,
            T,
            1.0,
        ]
        self.wt_gmix = ngmix.GMixModel(wtpars, "gauss")

    def get_moments(self, pars, more=False):
        """
        pars are [cen1, cen2, shear1, shear2]
        """
        obs=self.obs
        obs.jacold=obs.jacobian

        row_offset, col_offset = pars[0],pars[1]
        s1, s2 = pars[2], pars[3]
        try:
            
            wt_gmix=self.wt_gmix
            try:
                shear=ngmix.Shape(s1, s2)

                j=util.get_sheared_jacobian(
                    self.jacobian,
                    shear,
                )

                obs.jacobian=j

                wt_gmix.set_cen(row_offset, col_offset)

                res=wt_gmix.get_weighted_moments(obs)

                if more:
                    retval=res

                else:
                    # these are the un-normalized moments
                    retval = res['pars'][0:0+4]

            except ngmix.GMixRangeError as err:
                retval = zeros(4) + 9.9e9
        finally:
            obs.jacobian=obs.jacold

        return retval

class NullerGaussShape(NullerBase):
    def __init__(self, obs, T, cen, **kw):

        super(NullerGaussShape,self).__init__(obs, **kw)
        self.npars=2

        self.T = T

        wtpars = [
            cen[0],
            cen[1],
            0.0,
            0.0,
            T,
            1.0,
        ]
        self.wt_gmix = ngmix.GMixModel(wtpars, "gauss")

    def get_moments(self, pars, more=False):
        """
        pars are [cen1, cen2, shear1, shear2]
        """
        obs=self.obs
        obs.jacold=obs.jacobian

        s1, s2 = pars[0], pars[1]
        try:
            
            wt_gmix=self.wt_gmix
            try:
                shear=ngmix.Shape(s1, s2)

                j=util.get_sheared_jacobian(
                    self.jacobian,
                    shear,
                )

                obs.jacobian=j

                res=wt_gmix.get_weighted_moments(obs)

                if more:
                    retval=res

                else:
                    # these are the un-normalized moments
                    retval = res['pars'][2:2+2]

            except ngmix.GMixRangeError as err:
                retval = zeros(2) + 9.9e9
        finally:
            obs.jacobian=obs.jacold

        return retval


class MetacalNullGaussFitter(SimpleFitterBase):

    def _setup(self, *args, **kw):
        super(MetacalNullGaussFitter,self)._setup(*args, **kw)
        self['npars']=4

        self.metacal_types=self['metacal_pars'].get('types',DEFTYPES)

    def _dofit(self, imdict):

        obs=imdict['obs']
        odict=self._get_metacal(obs)


        guess0, T = self._get_weight_and_guess_from_admom(obs)

        res=self._do_metacal(odict, guess0, T)

        if self['metacal_pars']['do_noise_shear']:
            #self._fit_lines(res)
            nobs=ngmix.simobs.simulate_obs(None, obs)
            ndict=self._get_metacal(nobs)
            #self._add_noise_shears(ndict, T, res)

            nguess0 = guess0.copy()
            nguess0[2:2+2] = 0.0
            nres=self._do_metacal(ndict, nguess0, T, shape_only=True)
            for type in res:
                gn=nres[type]['g']
                res[type]['g_noise'] = gn
                if type=='1p':
                    print("g_noise:",gn[0], gn[1])

        res['flags']=0
        return res

    def _fit_lines(self, res):
        """
        we want to predict gamma, so the lines are

            gamma1 = A1 * M1 + B1
            gamma2 = A2 * M2 + B2
        
        """
        nuse=3

        g1off=array([-0.01, 0.0, 0.01])
        g2off=array([-0.01, 0.0, 0.01])

        M1 = zeros(nuse)
        M2 = zeros(nuse)

        for type in res:
            tres=res[type]

            tpars = tres['pars']

            nuller=tres['nuller']

            send_pars=tpars.copy()

            g1 = tpars[2] + g1off
            g2 = tpars[3] + g2off

            for i in xrange(nuse):

                send_pars[2] = g1[i]
                send_pars[3] = g2[i]

                if send_pars[2]**2 + send_pars[3]**2 > 0.999:
                    raise TryAgainError("could not fit line")

                tmoms = nuller.get_moments(send_pars)
                M1[i] = tmoms[2]
                M2[i] = tmoms[3]


            a1, b1 = util.fit_line(g1, M1)
            a2, b2 = util.fit_line(g2, M2)

            tres['coeffs1'] = (a1, b1)
            tres['coeffs2'] = (a2, b2)

            if False:
                import biggles
                tab=biggles.Table(1,2)
                plt1=biggles.plot(M1, g1, visible=False)
                plt1=biggles.plot(M1, a1*M1 + b1, type='solid', color='blue', plt=plt1, visible=False)

                plt2=biggles.plot(M2, g2, visible=False)
                plt2=biggles.plot(M2, a2*M2 + b2, type='solid', color='blue', plt=plt2, visible=False)

                tab[0,0] = plt1
                tab[0,1] = plt2
                tab.aspect_ratio=0.5

                tab.show(width=800,height=800)
                if 'q'==raw_input("hit a key: "):
                    stop

    def _add_noise_shears(self, ndict, T, res):

        maxiter=self['max_pars']['lm_pars']['maxfev']
        
        for type in res:
            nobs = ndict[type]

            tres=res[type]
            coeffs1 = tres['coeffs1']
            coeffs2 = tres['coeffs2']

            nuller=NullerGauss(nobs, T, maxiter=maxiter)

            moms=nuller.get_moments(tres['pars'])
            s1 = coeffs1[0]*moms[2] + coeffs1[1]
            s2 = coeffs2[0]*moms[3] + coeffs2[1]

            res[type]['g_noise'] = (s1,s2)

    def _get_metacal(self, obs):
        mcpars=self['metacal_pars']

        odict=ngmix.metacal.get_all_metacal(
            obs,
            rng=self.rng,
            **mcpars
        )

        return odict


    def _do_metacal(self, odict, guess0, T, shape_only=False):
        res={}

        for type in self.metacal_types:
            #print("    doing metacal:",type)
            obs=odict[type]

            res[type] = self._do_one_null_fit(obs, guess0, T, shape_only=shape_only)

        return res

    def _do_one_null_fit(self, obs, guess0, T, shape_only=False):

        mpars=self['max_pars']

        if shape_only:
            cen=guess0[0:0+2]
            guess0=guess0[2:2+2]
            nuller=NullerGaussShape(obs, T, cen, maxiter=mpars['lm_pars']['maxfev'])
        else:
            nuller=NullerGauss(obs, T, maxiter=mpars['lm_pars']['maxfev'])

        for i in xrange(mpars['ntry']):

            guess=self._get_guess(guess0)
            nuller.go(guess)

            res=nuller.get_result()
            if res['flags']==0:
                break

        if res['flags'] != 0:
            raise TryAgainError("        failed to find null")

        #print_pars(res['pars'], front="pars:   ")
        res['ntry'] = i+1

        res['nuller'] = nuller

        return res

    def _get_guess(self, guess0):

        guess=guess0.copy()

        if guess.size == 2:
            sh=self.get_shape_guess(guess[0], guess[1])
            guess[:] = (sh.g1,sh.g2)
        else:
            guess[0:0+2] += self.rng.uniform(low=-0.1,high=0.1,size=2)
            sh=self.get_shape_guess(guess[0], guess[1])
            guess[2:2+2] = (sh.g1, sh.g2)

        return guess

    def get_shape_guess(self, g1, g2, width=0.01, nmaxrand=10):
        """
        Get guess, making sure in range
        """
        rng=self.rng

        g=sqrt(g1**2 + g2**2)
        if g > 1.0:
            g1,g2=rng.uniform(low=-0.1, high=0.1, size=2)

        shape=ngmix.Shape(g1, g2)
        shape_guess=None

        for i in xrange(nmaxrand):
            try:
                g1_offset = rng.uniform(low=-width,high=width)
                g2_offset = rng.uniform(low=-width,high=width)
                shape_new=shape.get_sheared(g1_offset, g2_offset)

                if shape_new.g < 0.98:
                    shape_guess=shape_new
                    break

            except GMixRangeError:
                pass

        if shape_guess is None:
            g1g=rng.uniform(low=-0.1, high=0.1)
            g2g=rng.uniform(low=-0.1, high=0.1)
            shape_guess=ngmix.Shape(g1g,g2g)

        return shape_guess

    def _get_weight_and_guess_from_admom(self, obs):
        """
        measure adaptive moments
        """

        try:
            ampars=self['admom_pars']
            ntry=ampars.pop('ntry',4)

            fitter=ngmix.admom.Admom(obs, **ampars)

            for i in xrange(ntry):
                guess=self._get_admom_guess(obs)

                fitter.go(guess)

                res=fitter.get_result()

                if res['flags'] == 0:
                    break

            if res['flags'] != 0:
                raise TryAgainError("        initial admom failed")

            gmix=fitter.get_gmix()
            c1,c2=gmix.get_cen()
            g1,g2,T0 = gmix.get_g1g2T()

            guess = numpy.array( [c1,c2,-g1,-g2] )
            #T = ngmix.moments.get_Tround(T, g1, g2)
            T = T0

        except ngmix.GMixRangeError:
            raise TryAgainError("        failed to get guess from admom")

        return guess, T

    def _get_admom_guess(self, obs):

        rng=self.rng

        scale=obs.jacobian.get_scale()
        Tguess = 10.0*scale**2

        pars=zeros(6)
        pars[0:0+2] = rng.uniform(low=-0.5*scale, high=0.5*scale, size=2)
        pars[2:2+2] = rng.uniform(low=-0.3, high=0.3, size=2)
        pars[4]     = Tguess*(1.0 + rng.uniform(low=-0.1, high=0.1))
        pars[5]     = 1.0

        guess=ngmix.GMixModel(pars, "gauss")
        return guess



    def _fit_gauss(self, obs):
        """
        fit a gaussian to get a good weight function
        """

        mpars=self['max_pars']

        Tguess=4.0
        runner = ngmix.bootstrap.PSFRunner(
            obs,
            "gauss",
            Tguess,
            mpars['lm_pars'],
        )
        runner.go(ntry=mpars['ntry'])

        res=runner.fitter.get_result()
        if res['flags'] != 0:
            raise TryAgainError("failed to fit gaussian")

        pars=res['pars']

        return pars


    def _print_res(self,allres):
        """
        print some stats
        """

        res=allres['noshear']
        if 'nfev' in res:
            mess="    s2n: %.1f  ntry: %d  nfev: %d"
            mess = mess % (res['s2n'],res['ntry'],res['nfev'])
            print(mess)

        print_pars(res['pars'],      front='        pars: ')

        if allres['pars_true'][0] is not None:
            print_pars(allres['pars_true'], front='        true: ')

    def _get_dtype(self):
        """
        get the dtype for the output struct
        """
        npars=self['npars']

        # super of super
        dt=super(MetacalNullGaussFitter,self)._get_dtype()


        do_noise_shear=self['metacal_pars']['do_noise_shear']

        for type in self.metacal_types:

            if type=='noshear':
                back=''
            else:
                back='_%s' % type

            dt += [
                ('mcal_pars%s' % back,'f8',npars),

                ('mcal_g%s' % back,'f8',2),
                ('mcal_g_cov%s' % back,'f8',(2,2)),

                ('mcal_s2n%s' % back,'f8'),

                ('mcal_sums%s' % back,'f8',6),

                ('mcal_nfev%s' % back, 'i4'),
                ('mcal_ntry%s' % back, 'i4'),
            ]

            if do_noise_shear and type != 'noshear':
                dt += [('mcal_g_noise%s' % back, 'f8',2)]

            if type=='noshear':
                dt += [
                    ('mcal_wsum','f8'),
                ]

            dt += [
            ]

        return dt


    def _copy_to_output(self, res, i):
        """
        copy parameters specific to this class
        """

        # note copying super of our super, since
        # we didn't do a regular fit
        super(SimpleFitterBase,self)._copy_to_output(res, i)

        do_noise_shear=self['metacal_pars']['do_noise_shear']

        d=self.data

        for type in self.metacal_types:

            tres=res[type]
            if type=='noshear':
                back=''
            else:
                back='_%s' % type

            d['mcal_pars%s' % back][i] = tres['pars']

            d['mcal_g%s' % back][i] = tres['g']
            d['mcal_g_cov%s' % back][i] = tres['g_cov']

            d['mcal_s2n%s' % back][i] = tres['s2n']

            d['mcal_sums%s' % back][i] = tres['sums']

            d['mcal_nfev%s' % back][i] = tres['nfev']
            d['mcal_ntry%s' % back][i] = tres['ntry']

            if do_noise_shear and type != 'noshear':
                d['mcal_g_noise%s' % back][i] = tres['g_noise']

            if type=='noshear':
                for p in ['wsum']:

                    if p in tres:
                        name='mcal_%s' % p
                        d[name][i] = tres[p]


def _get_weight(obs, config):
    """
    for now used fixed gaussian weight function
    """
    import galsim
    from math import sqrt

    wtpars=config['weight']
    if 'r50' in wtpars:
        #sigma = wtpars['r50']/sqrt(2(log(2)))
        sigma = wtpars['r50']/1.1774100225154747
        T=2*sigma**2

    elif 'sigma' in wtpars:
        sigma = wtpars['sigma']
        T=2*sigma**2

    elif 'T' in wtpars:
        T = wtpars['T']

    else:
        raise RuntimeError("have either r50 or sigma "
                           "in weight section")

    # for now at the jacobian center
    pars=[0.0, 0.0, 0.0, 0.0, T, 1.0]
    wt_gmix=ngmix.GMixModel(pars, "gauss")

    return wt_gmix

def _M1M2T_to_matrix(M1, M2, T):
    Irr, Irc, Icc = _M1M2T_to_IrrIrcIcc(M1, M2, T)

    M = zeros( (2,2) )
    M[0,0] = Irr
    M[0,1] = Irc
    M[1,0] = Irc
    M[1,1] = Icc
    return M


def _M1M2T_to_IrrIrcIcc(M1, M2, T):
    Icc = 0.5*(T + M1)
    Irr = 0.5*(T - M1)
    Irc = 0.5*M2

    return Irr, Irc, Icc


def _doadd_single_obs_combined(obs, nobs, nrand):
    im  = obs.image
    nim = nobs.image

    image = im + nim

    weight=numpy.zeros(obs.weight.shape)

    wpos=numpy.where(
        (obs.weight != 0.0) &
        (nobs.weight != 0.0)
    )
    if wpos[0].size > 0:
        tvar = obs.weight*0
        # add the variances
        tvar[wpos] = (
            nrand/obs.weight[wpos]  +
            1.0/nobs.weight[wpos]
        )
        weight[wpos] = 1.0/tvar[wpos]

    newobs=ngmix.Observation(
        image,
        weight=weight,
        jacobian=obs.jacobian.copy(),
    )
    return newobs

def get_all_metacal_multi_combined(obs, step=0.01, **kw):
    """
    do multiple random realizations
    currently only for single obs input
    """
    assert isinstance(obs, ngmix.Observation)

    nrand=kw.pop('nrand',1)
    print("doing nrand:",nrand)

    orig_obsdict = ngmix.metacal._get_all_metacal(obs, step=step, **kw)

    obsdict={}
    for key in orig_obsdict:
        obsdict[key] = ngmix.ObsList()
    for i in xrange(nrand):
        # Using None for the model means we get just noise
        noise_obs = ngmix.simobs.simulate_obs(None, obs, **kw)

        # rotate by 90
        ngmix.metacal._rotate_obs_image(noise_obs, k=1)

        tnoise_obsdict = ngmix.metacal._get_all_metacal(noise_obs, step=step, **kw)

        for key in tnoise_obsdict:
            tobs=orig_obsdict[key]

            nobs=tnoise_obsdict[key]
            ngmix.metacal._rotate_obs_image(nobs, k=3)

            newobs = _doadd_single_obs_combined(tobs, nobs, nrand)

            obsdict[key].append( newobs )

    return obsdict

def _doadd_single_obs(obs, nobs):
    im  = obs.image
    nim = nobs.image

    image = im + nim

    weight=numpy.zeros(obs.weight.shape)

    wpos=numpy.where(
        (obs.weight != 0.0) &
        (nobs.weight != 0.0)
    )
    if wpos[0].size > 0:
        tvar = obs.weight*0
        # add the variances
        tvar[wpos] = (
            1.0/obs.weight[wpos]  +
            1.0/nobs.weight[wpos]
        )
        weight[wpos] = 1.0/tvar[wpos]

    newobs=ngmix.Observation(
        image,
        weight=weight,
        jacobian=obs.jacobian.copy(),
    )
    return newobs


def get_all_metacal_multi(obs, step=0.01, **kw):
    """
    do multiple random realizations
    currently only for single obs input
    """
    assert isinstance(obs, ngmix.Observation)

    nrand=kw.pop('nrand',1)
    print("doing nrand:",nrand)

    orig_obsdict = ngmix.metacal._get_all_metacal(obs, step=step, **kw)

    obsdict_list=[]

    for i in xrange(nrand):
        # Using None for the model means we get just noise
        noise_obs = ngmix.simobs.simulate_obs(None, obs, **kw)

        # rotate by 90
        ngmix.metacal._rotate_obs_image(noise_obs, k=1)

        noise_obsdict = ngmix.metacal._get_all_metacal(noise_obs, step=step, **kw)

        obsdict={}
        for key in noise_obsdict:
            tobs=orig_obsdict[key]

            nobs=noise_obsdict[key]
            ngmix.metacal._rotate_obs_image(nobs, k=3)

            newobs = _doadd_single_obs(tobs, nobs)

            obsdict[key] = newobs

        obsdict_list.append(obsdict)

    return obsdict_list

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
