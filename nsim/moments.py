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

        return res

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
                    ('mcal_gpsf','f8',2),
                    ('mcal_Tpsf','f8'),
                ]

            dt += [
                #('mcal_s2n_r%s' % back,'f8'),
                #('mcal_s2n_w%s' % back,'f8'),

                #('mcal_r50%s' % back,'f8'),
                #('mcal_r50_s2n%s' % back,'f8'),
                ('mcal_s2n%s' % back,'f8'),
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

            d['mcal_s2n%s' % back][i] = tres['s2n']

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
        print("    e1e2:  %g %g" % tuple(subres['g']))
        print("    e_err: %g" % numpy.sqrt(subres['g_cov'][0,0]))

        print_pars(subres['pars'],      front='        pars: ')

        print('        true: ', res['pars_true'])



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

    def _get_guess(self):
        rng=self.rng

        pars=zeros(6)
        pars[0:0+2] = rng.uniform(low=-0.5, high=0.5, size=2)
        pars[2:2+2] = rng.uniform(low=-0.3, high=0.3, size=2)
        pars[4]     = 10.0*(1.0 + rng.uniform(low=-0.1, high=0.1))
        pars[5]     = 1.0

        guess=ngmix.GMixModel(pars, "gauss")
        return guess

class MetacalMomentsAM(MetacalMomentsFixed):
    def _dofit(self, imdict):

        obs=imdict['obs']

        self._set_metacal(obs)

        res=self._do_metacal()
        res['flags']=0
        return res

    def _measure_moments(self, obs):
        """
        moments of the power spectrum

        todo: subtract noise ps. Deal with errors
        """

        ampars=self['admom_pars']
        ntry=ampars.pop('ntry',4)

        fitter=ngmix.admom.Admom(obs, **ampars)

        for i in xrange(ntry):
            guess=self._get_guess()

            fitter.go(guess)
            res=fitter.get_result()

            if res['flags'] == 0:
                break

        if res['flags'] != 0:
            raise TryAgainError("admom failed")

        res['g']     = res['e']
        res['g_cov'] = numpy.diag( [ res['e_err']**2 ]*2 )
        res['pars']  = res['sums']

        # not right pars cov
        res['pars_cov']=res['sums_cov']*0 + 9999.e9

        res['s2n']      = res['s2n']
        res['flux']     = res['flux']
        res['flux_s2n'] = res['s2n']

        return res

    def _copy_to_output(self, res, i):
        """
        copy parameters specific to this class
        """

        # note copying super of our super, since
        # we didn't do a regular fit
        super(MetacalMomentsAM,self)._copy_to_output(res, i)

        d=self.data

        for type in self.metacal_types:
            if type=='noshear':
                back=''
            else:
                back='_%s' % type

            tres=res[type]
            d['mcal_numiter%s' % back][i] = tres['numiter']



    def _get_dtype(self):
        """
        get the dtype for the output struct
        """
        # super of super
        dt=super(MetacalMomentsAM,self)._get_dtype()

        for type in self.metacal_types:
            if type=='noshear':
                back=''
            else:
                back='_%s' % type

            dt += [('mcal_numiter%s' % back,'i4')]

        return dt

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
            guess=self._get_guess()

            fitter.go(guess)
            res=fitter.get_result()

            if res['flags'] == 0:
                break

        if res['flags'] != 0:
            raise TryAgainError("admom failed for psf")

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
            raise TryAgainError("bad moments")

        wpars=res['pars']
        M1 = wpars[2]/wpars[5]
        M2 = wpars[3]/wpars[5]
        T  = wpars[4]/wpars[5]

        Mmeas = _M1M2T_to_matrix(M1, M2, T)

        try:
            Mmeas_inv = numpy.linalg.pinv(Mmeas)
        except numpy.linalg.LinAlgError:
            raise TryAgainError("could not invert observed moment matrix")


        Minv_deweight = Mmeas_inv - self.Mwt_inv

        try:
            Mdeweight = numpy.linalg.pinv(Minv_deweight )
        except numpy.linalg.LinAlgError:
            raise TryAgainError("could not invert deweighted moment matrix")
        
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

        return res


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

