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

        mpars=self['metacal_pars']
        deftypes=[
            'noshear',
            '1p','1m','2p','2m',
        ]
        sym=mpars.get('symmetrize_psf',False)
        if not sym:
            deftypes += [
                '1p_psf','1m_psf',
                '2p_psf','2m_psf',
            ]
        self.metacal_types=mpars.get('types',deftypes)
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

    def _get_metacal(self, obs):
        mcpars=self['metacal_pars']
        analytic_psf=None
        odict=ngmix.metacal.get_all_metacal(
            obs,
            psf=analytic_psf,
            rng=self.rng,
            **mcpars
        )

        return odict

    def _do_metacal(self, odict):
        res={}

        for type,obs in odict.iteritems():

            tres,fitter=self._measure_moments(obs)
            if tres['flags'] != 0:
                raise TryAgainError("        bad T")


            if type=='noshear':
                #pres,fitter=self._measure_moments(obs.psf_nopix)
                pres,fitter=self._measure_moments(obs.psf)
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

            if type=='noshear':
                for p in ['pars_cov','wsum','psfrec_g','psfrec_T']:

                    if p in tres:
                        name='mcal_%s' % p
                        d[name][i] = tres[p]

    def _print_res(self,res):
        """
        print some stats
        """

        preres=res['prefit']
        subres=res['noshear']

        print()
        print("    s2n: %g" % preres['s2n'])

        s2n=subres['s2n']
        g=subres['g']
        gerr=diag(sqrt(subres['g_cov']))

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

class MetacalMomentsAM(MetacalMomentsFixed):
    def _dofit(self, imdict):


        obs=imdict['obs']
        if False:
            self._do_plots(obs)

        psfres, fitter = self._measure_moments(obs.psf)

        obsdict=self._get_metacal(obs)

        if False:
            tobs = obsdict['noshear']
            noshear_obs = ngmix.Observation(
                tobs.image_orig,
                tobs.weight_orig,
                jacobian=tobs.jacobian.copy(),
            )

            pre_res,self.pre_fitter=self._measure_moments(noshear_obs)
        else:
            pre_res,self.pre_fitter=self._measure_moments(obs)


        pre_res['psfrec_g']= psfres['g']
        pre_res['psfrec_T']= psfres['T']


        res=self._do_metacal(obsdict)
        res['prefit'] = pre_res

        res['flags']=0
        return res

    def _measure_moments(self, obs):
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

            if res['flags'] == 0:
                break

        if res['flags'] != 0:
            raise TryAgainError("        admom failed")

        self._set_some_pars(res)
        return res, fitter

    def _set_some_pars(self, res):
        res['g']     = res['e']
        res['g_cov'] = res['e_cov']

        # not right pars cov
        res['pars_cov']=res['sums_cov']*0 + 9999.e9
        res['flux_s2n'] = res['s2n']


    def _copy_to_output(self, res, i):
        """
        copy parameters specific to this class
        """

        # note copying super of our super, since
        # we didn't do a regular fit
        super(MetacalMomentsAM,self)._copy_to_output(res, i)

        d=self.data

        pres=res['prefit']
        d['pars'][i] = pres['pars']
        d['g'][i] = pres['g']
        d['g_cov'][i] = pres['g_cov']
        d['s2n'][i] = pres['s2n']
    
        d['psfrec_g'][i] = pres['psfrec_g']
        d['psfrec_T'][i] = pres['psfrec_T']

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

        dt += [
            ('psfrec_g','f8',2),
            ('psfrec_T','f8'),
        ]
        for type in self.metacal_types:
            if type=='noshear':
                back=''
            else:
                back='_%s' % type

            dt += [('mcal_numiter%s' % back,'i4')]

        return dt


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

class AMFitter(MetacalMomentsAM):
    def _setup(self, *args, **kw):
        super(MetacalMomentsFixed,self)._setup(*args, **kw)

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
        analytic_psf=None
        odict_list=get_all_metacal_multi(
            obs,
            psf=analytic_psf,
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

