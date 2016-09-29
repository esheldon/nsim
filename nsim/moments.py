"""
moment based code
"""
from __future__ import print_function

import numpy
from numpy import array, zeros, ones, log, log10, exp, sqrt, diag
from numpy import where, isfinite

import ngmix
from ngmix.fitting import print_pars
from .util import TryAgainError

from .fitters import SimpleFitterBase

class MetacalMomentsFixed(SimpleFitterBase):
    """
    - might be effects not captured by the galsim
      Convolve
    - probably should dilate while we are symmetrizing
    """
    def _setup(self, *args, **kw):
        super(MetacalMomentsFixed,self)._setup(*args, **kw)

        self.metacal_types=[
            'noshear',
            '1p','1m','2p','2m',
        ]

    def _dofit(self, imdict):

        obs=imdict['obs']

        self._set_weight(obs) 

        self._set_metacal(obs)

        #self._set_F(self.kmcal.dim)

        res=self._do_metacal()
        res['flags']=0
        return res

    def _set_metacal(self, obs):
        mcpars=self['metacal_pars']
        analytic_psf=None
        self.odict=ngmix.metacal.get_all_metacal(
            obs,
            psf=analytic_psf,
            rng=self.rng,
            **mcpars
        )

    def _do_metacal(self):
        res={}

        for type,obs in self.odict.iteritems():

            tres=self._measure_moments(obs)
            if tres['flags'] != 0:
                raise TryAgainError("bad T")

            res[type]=tres

        return res

    def _measure_moments(self, obs):
        """
        moments of the power spectrum

        todo: subtract noise ps. Deal with errors
        """

        jac=obs.jacobian
        cen=jac.get_cen()

        scale=jac.get_scale()
        rmax=min(cen)*scale

        # flags would only be set if there were ivar <= 0
        res=self.wt_gmix.get_weighted_moments(obs, rmax=rmax)

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

        return res

    def _set_weight(self, obs):
        """
        for now used fixed gaussian weight function
        """
        import galsim
        from math import sqrt

        wtpars=self['weight']
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
        self.wt_gmix=ngmix.GMixModel(pars, "gauss")

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
                    ('mcal_gpsf','f8',2),
                    ('mcal_Tpsf','f8'),
                ]

            dt += [
                #('mcal_s2n_r%s' % back,'f8'),
                #('mcal_s2n_w%s' % back,'f8'),

                #('mcal_r50%s' % back,'f8'),
                #('mcal_r50_s2n%s' % back,'f8'),
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

            #d['mcal_s2n_w%s' % back][i] = tres['s2n_w']
            #d['mcal_s2n_r%s' % back][i] = tres['s2n_r']

            #r50 = tres['pars'][4]
            #r50_s2n = r50/sqrt(tres['pars_cov'][4,4])
            #d['mcal_r50%s' % back][i] = r50
            #d['mcal_r50_s2n%s' % back][i] = r50_s2n

            d['mcal_flux%s' % back][i] = tres['flux']
            d['mcal_flux_s2n%s' % back][i] = tres['flux_s2n']

            if type=='noshear':
                for p in ['pars_cov','wsum','gpsf','Tpsf']:

                    if p in tres:
                        name='mcal_%s' % p
                        d[name][i] = tres[p]

    def _print_res(self,res):
        """
        print some stats
        """

        subres=res['noshear']

        print("    flux s2n: %g" % subres['flux_s2n'])
        print("    g1g2: %g %g" % tuple(subres['g']))

        print_pars(subres['pars'],      front='        pars: ')
        print_pars(subres['pars_err'],  front='        perr: ')

        print_pars(res['pars_true'], front='        true: ')



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


