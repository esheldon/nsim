from __future__ import print_function
import sys
import os
import shutil
from pprint import pprint
import numpy
from numpy import sqrt, array, diag, ones, zeros
from numpy import where, newaxis, exp, log
from numpy import newaxis
import yaml

import fitsio

import ngmix
import nsim
from . import files
from . import shearpdf
from .util import Namer

import argparse
import esutil as eu
from esutil.numpy_util import between

try:
    import reredux
except ImportError:
    pass

class Summer(dict):
    def __init__(self, conf, args, shears=None, nshear=None):
        self.update(conf)
        self.args=args
        self.chunksize=args.chunksize

        self._set_select()

        self.namer=Namer(front='mcal')
        self.gpsf_name='mcal_gpsf'

        self.step = self['metacal_pars'].get('step',0.01)

        self.shears = shears

        if shears is None and nshear is None:
            raise RuntimeError("send either shears= or nshear=")

        if shears is not None:
            self['nshear']=len(self.shears)
        else:
            self['nshear'] = nshear

    def go(self):

        if self.args.fit_only:
            self.means=self._read_means()
        else:

            extra=self._get_fname_extra()

            sums = self.do_sums()

            args=self.args

            g,gpsf,R,Rpsf,Rsel,Rsel_psf=self._average_sums(sums)

            means=get_mean_struct(self['nshear'])
            means_nocorr=get_mean_struct(self['nshear'])
            
            wkeep,=numpy.where(sums['wsum'] > 0)

            for i in xrange(self['nshear']):

                if self.shears is not None:
                    shear_true = self.shears[i]
                else:
                    shear_true = sums['shear_true'][i]

                gmean = g[i]

                c        = (Rpsf+Rsel_psf)*gpsf[i]
                c_nocorr = Rpsf*gpsf[i]

                shear        = (gmean-c)/(R+Rsel)
                shear_nocorr = (gmean-c_nocorr)/R

                means['shear'][i] = shear
                means['shear_err'][i] = 1.0
                if isinstance(shear_true,ngmix.Shape):
                    means['shear_true'][i,0] = shear_true.g1
                    means['shear_true'][i,1] = shear_true.g2
                else:
                    means['shear_true'][i] = shear_true

                means_nocorr['shear'][i] = shear_nocorr
                means_nocorr['shear_err'][i] = 1.0
                means_nocorr['shear_true'][i] = means['shear_true'][i]

            means=means[wkeep]
            means_nocorr=means_nocorr[wkeep]
            self.means=means
            self.means_nocorr=means_nocorr
            self._write_means()

        if self.do_selection:
            print("without correction")
            junk=reredux.averaging.fit_m_c(self.means_nocorr)
            junk=reredux.averaging.fit_m_c(self.means_nocorr,onem=True)
            print("\nwith correction")

        self.fits=reredux.averaging.fit_m_c(self.means)
        self.fitsone=reredux.averaging.fit_m_c(self.means,onem=True)


    def get_run_output(self, run):
        """
        collated file
        """

        fname = nsim.files.get_output_url(run, 0, 0)
        if self.args.cache:
            origname=fname
            bname = os.path.basename(origname)
            fname=os.path.join('$TMPDIR', bname)
            fname=os.path.expandvars(fname)
            if not os.path.exists(fname):
                print("copying to cache: %s -> %s" % (origname,fname))
                shutil.copy(origname, fname)

        return fname

    def do_sums(self):

        args=self.args
        chunksize=self.chunksize
        n=self.namer

        sums=None
        ntot=0
        for run in args.runs:
            if args.ntest is not None and ntot > args.ntest:
                break

            fname=self.get_run_output(run)
            print(fname)
            with fitsio.FITS(fname) as fits:

                hdu=fits[1]

                nrows=hdu.get_nrows()
                nchunks = nrows//chunksize

                if (nrows % chunksize) > 0:
                    nchunks += 1

                beg=0
                for i in xrange(nchunks):
                    print("    chunk %d/%d" % (i+1,nchunks))

                    end=beg+chunksize

                    data = hdu[beg:end]

                    data=self._preselect(data)

                    ntot += data.size

                    #sums=self.do_sums1(data, sums=sums)
                    if False and n('g') not in data.dtype.names:
                        if True:
                            sums=self.do_sums1_moms_wt(data, sums=sums)
                        elif False:
                            sums=self.do_sums1_moms(data, sums=sums)
                        elif False:
                            sums=self.do_sums1_moms_psfcorr(data, sums=sums)
                    else:
                        sums=self.do_sums1(data, sums=sums)

                    beg = beg + chunksize

                    if args.ntest is not None and ntot > args.ntest:
                        break

        return sums

    def _preselect(self, data):
        """
        sub-classes might make a pre-selection, e.g. of some flags
        """

        if False:
            g_1p = data['mcal_g_1p'][:,0]
            g_1m = data['mcal_g_1m'][:,0]
            g_2p = data['mcal_g_2p'][:,1]
            g_2m = data['mcal_g_2m'][:,1]

            R1=(g_1p-g_1m)/(2.0*self.step)
            R2=(g_2p-g_2m)/(2.0*self.step)

            w,=numpy.where(between(R1, -5, 6.5) & between(R2, -5, 6.5))
            print("kept %d/%d in preselect" % (w.size, data.size))
            data=data[w]

        return data


    def _get_weights(self, data, w, type):

        n=self.namer
        if self.args.weighted:
            if type=='noshear':
                name=n('g_cov')
            else:
                name=n('g_cov_%s' % type)

            g_cov=data[name][w]

            wts=get_noise_weights(g_cov, self.args)

        else:
            wts=numpy.ones(w.size)

        wa = wts[:,numpy.newaxis]
        return wts, wa

    def _get_g_name(self, data, type):
        n=self.namer
        if type=='noshear':
            name=n('g')
        else:
            name=n('g_%s' % type)

        return name

    def _get_g(self, data, w, type):
        name=self._get_g_name(data, type)

        if name not in data.dtype.names:
            g = None
        else:
            g = data[name][w]

        return g

    def _get_gpsf(self, data, w):
        return data[self.gpsf_name][w]

    def do_sums1(self, data, sums=None):
        """
        just a binner and summer, no logic here
        """


        nshear=self['nshear']
        args=self.args

        h,rev = eu.stat.histogram(data['shear_index'],
                                  min=0,
                                  max=nshear-1,
                                  rev=True)
        nind = h.size
        assert nshear==nind

        if sums is None:
            sums=self._get_sums_struct()

        ntot=0
        nkeep=0
        wtmax=0.0
        wttot=0.0
        for i in xrange(nshear):
            if rev[i] != rev[i+1]:
                wfield=rev[ rev[i]:rev[i+1] ]

                # first select on the noshear measurement
                if self.select is not None:
                    w=self._do_select(data, wfield)
                    w=wfield[w]
                else:
                    w=wfield

                if w.size == 0:
                    continue

                ntot  += wfield.size
                nkeep += w.size

                if 'shear_true' in data.dtype.names:
                    # should all be the same, so just copy
                    # the first one.  We will end up copying over this
                    # each time, but that's ok
                    sums['shear_true'][i] += data['shear_true'][w[0]]

                g = self._get_g(data, w, 'noshear')
                wts, wa = self._get_weights(data, w, 'noshear')

                wttot += wts.sum()
                twtmax = wts.max()
                if twtmax > wtmax:
                    wtmax = twtmax

                sums['g'][i]    += (g*wa).sum(axis=0)
                gpsf=self._get_gpsf(data, w)
                sums['gpsf'][i] += (gpsf*wa).sum(axis=0)
                sums['wsum'][i] += wts.sum()

                for type in ngmix.metacal.METACAL_TYPES:
                    if type=='noshear':
                        continue

                    sumname='g_%s' % type

                    g=self._get_g(data, w, type)

                    if g is not None:
                        # using the same weights, based on unsheared
                        # parameters
                        sums[sumname][i] += (g*wa).sum(axis=0)

                # now the selection terms
                if self.do_selection is not None:
                    for type in ngmix.metacal.METACAL_TYPES:
                        if type=='noshear':
                            continue

                        g_name=self._get_g_name(data, type)

                        if g_name in data.dtype.names:

                            wsumname = 's_wsum_%s' % type
                            sumname = 's_g_%s' % type

                            if self.select is not None:
                                w=self._do_select(data, wfield, type)
                                w=wfield[w]
                            else:
                                w=wfield

                            # weights based on sheared parameters
                            g=self._get_g(data, w, 'noshear')
                            wts,wa=self._get_weights(data, w, type)

                            sums[sumname][i] += (g*wa).sum(axis=0)
                            sums[wsumname][i] += wts.sum()

        if self.args.weighted:
            effnum=wttot/wtmax/ntot
            print("        effective fraction: %g" % effnum)

        if self.select is not None:
            self._print_frac(ntot,nkeep)
        return sums

    def _get_M1M2corr(self, pars, psf_pars):
        M1 = pars[:,2]
        M2 = pars[:,3]
        T = pars[:,4]
        Fsum = pars[:,5]

        psfM1 = psf_pars[:,2]
        psfM2 = psf_pars[:,3]
        psfT = psf_pars[:,4]
        psfFsum = psf_pars[:,5]

        Frat = Fsum/psfFsum
        M1 = M1 - psfM1*Frat
        M2 = M2 - psfM2*Frat
        T = T - psfT*Frat

        w,=numpy.where(T == 0)
        if w.size > 0:
            print("found",w.size,"zeros")
        e1=M1/T
        e2=M2/T
        return e1,e2

    def _get_M1M2corr_nodiv(self, pars, psf_pars):
        M1 = pars[:,2]
        M2 = pars[:,3]
        T = pars[:,4]
        Fsum = pars[:,5]

        psfM1 = psf_pars[:,2]
        psfM2 = psf_pars[:,3]
        psfT = psf_pars[:,4]
        psfFsum = psf_pars[:,5]

        Frat = Fsum/psfFsum
        M1 = M1 - psfM1*Frat
        M2 = M2 - psfM2*Frat
        T = T - psfT*Frat

        return M1, M2




    def do_sums1_moms_psfcorr(self, data, sums=None):
        """
        just a binner and summer, no logic here
        """

        #getfunc=self._get_M1M2corr
        getfunc=self._get_M1M2corr_nodiv

        nshear=self['nshear']
        args=self.args

        h,rev = eu.stat.histogram(data['shear_index'],
                                  min=0,
                                  max=nshear-1,
                                  rev=True)
        nind = h.size
        assert nshear==nind

        if sums is None:
            sums=self._get_sums_struct()

        ntot=0
        nkeep=0
        for i in xrange(nshear):
            if rev[i] != rev[i+1]:
                wfield=rev[ rev[i]:rev[i+1] ]

                # first select on the noshear measurement
                if self.select is not None:
                    w=self._do_select(data['mcal_s2n'][wfield])
                    w=wfield[w]
                else:
                    w=wfield

                ntot  += wfield.size
                nkeep += w.size

                sums['wsum'][i] += w.size

                M1, M2 = getfunc(
                    data['mcal_pars'][w],
                    data['mcal_psf_pars'][w],
                )

                sums['g'][i,0] += M1.sum()
                sums['g'][i,1] += M2.sum()

                if 'mcal_gpsf' in data.dtype.names:
                    sums['gpsf'][i] += data['mcal_gpsf'][w].sum(axis=0)

                for type in ngmix.metacal.METACAL_TYPES:
                    if type=='noshear':
                        continue
                    mcalname='mcal_pars_%s' % type
                    psf_mcalname='mcal_psf_pars_%s' % type

                    if mcalname in data.dtype.names:
                        sumname='g_%s' % type

                        M1, M2 = getfunc(data[mcalname][w],data[psf_mcalname][w])
                        sums[sumname][i,0] += M1.sum()
                        sums[sumname][i,1] += M2.sum()
                    else:
                        #print("    skipping:",mcalname)
                        pass

                # now the selection terms

                if self.select is not None:
                    for type in ngmix.metacal.METACAL_TYPES:
                        if type=='noshear':
                            continue
                        s2n_name='mcal_s2n_%s' % type

                        if s2n_name in data.dtype.names:
                            wsumname = 's_wsum_%s' % type
                            sumname = 's_g_%s' % type

                            w=self._do_select(data[s2n_name][wfield])
                            w=wfield[w]
                            sums[wsumname][i] += w.size

                            M1, M2 = getfunc(data['mcal_pars'][w],data['mcal_psf_pars'][w])

                            sums[sumname][i,0]  += M1.sum()
                            sums[sumname][i,1]  += M2.sum()
                        else:
                            #print("    skipping:",s2n_name)
                            pass

        if self.select is not None:
            self._print_frac(ntot,nkeep)
        return sums


    def do_sums1_moms(self, data, sums=None):
        """
        just a binner and summer, no logic here
        """

        di=4

        s2n_name=self._get_s2n_name(data)

        nshear=self['nshear']
        args=self.args

        h,rev = eu.stat.histogram(data['shear_index'],
                                  min=0,
                                  max=nshear-1,
                                  rev=True)
        nind = h.size
        assert nshear==nind

        if sums is None:
            sums=self._get_sums_struct()

        ntot=0
        nkeep=0
        for i in xrange(nshear):
            if rev[i] != rev[i+1]:
                wfield=rev[ rev[i]:rev[i+1] ]

                # first select on the noshear measurement
                if self.select is not None:
                    w=self._do_select(data[s2n_name][wfield])
                    w=wfield[w]
                else:
                    w=wfield

                ntot  += wfield.size
                nkeep += w.size

                sums['wsum'][i] += w.size

                e=data['mcal_pars'][w,2:2+2]
                e[:,0] /= data['mcal_pars'][w,di]
                e[:,1] /= data['mcal_pars'][w,di]

                sums['g'][i] += e.sum(axis=0)

                if 'mcal_gpsf' in data.dtype.names:
                    sums['gpsf'][i] += data['mcal_gpsf'][w].sum(axis=0)

                for type in ngmix.metacal.METACAL_TYPES:
                    if type=='noshear':
                        continue
                    mcalname='mcal_pars_%s' % type

                    if mcalname in data.dtype.names:
                        sumname='g_%s' % type

                        e = data[mcalname][w,2:2+2]
                        e[:,0] /= data[mcalname][w,di]
                        e[:,1] /= data[mcalname][w,di]
                        sums[sumname][i] += e.sum(axis=0)
                    else:
                        #print("    skipping:",mcalname)
                        pass

                # now the selection terms

                if self.select is not None:
                    for type in ngmix.metacal.METACAL_TYPES:
                        if type=='noshear':
                            continue
                        ts2n_name='%s_%s' % (s2n_name,type)

                        if ts2n_name in data.dtype.names:
                            wsumname = 's_wsum_%s' % type
                            sumname = 's_g_%s' % type

                            w=self._do_select(data[ts2n_name][wfield])
                            w=wfield[w]
                            sums[wsumname][i] += w.size

                            e=data['mcal_pars'][w,2:2+2]
                            e[:,0] /= data['mcal_pars'][w,di]
                            e[:,1] /= data['mcal_pars'][w,di]
                            sums[sumname][i]  += e.sum(axis=0)
                        else:
                            #print("    skipping:",ts2n_name)
                            pass

        if self.select is not None:
            self._print_frac(ntot,nkeep)
        return sums

    def _get_noise_weight_mom(self, data, w, type):
        if type=='noshear':
            tstr=''
        else:
            tstr='_%s' % type
        parname='mcal_pars%s' % tstr
        covname='mcal_pars_cov%s' % tstr

        M1=data[parname][w,2]
        M2=data[parname][w,3]
        T=data[parname][w,4]

        e1=M1/T
        e2=M2/T

        VM1=data[covname][w,2,2]
        VM2=data[covname][w,3,3]
        VT=data[covname][w,4,4]
        Ve = VM1/T + VM2/T + VT*M1/T**2 + VT*M2/T**2

        Ti2=1/T**2
        Ve = ( VM1*Ti2 + e1**2 *VT * Ti2  +
               VM2*Ti2 + e2**2 *VT * Ti2 )
        wt = 1.0/(2*0.3**2 + Ve)

        return wt


    def do_sums1_moms_wt(self, data, sums=None):
        """
        just a binner and summer, no logic here
        """

        self.do_selection=True

        nshear=self['nshear']
        args=self.args

        h,rev = eu.stat.histogram(data['shear_index'],
                                  min=0,
                                  max=nshear-1,
                                  rev=True)
        nind = h.size
        assert nshear==nind

        if sums is None:
            sums=self._get_sums_struct()

        ntot=0
        nkeep=0
        for i in xrange(nshear):
            if rev[i] != rev[i+1]:
                wfield=rev[ rev[i]:rev[i+1] ]

                # first select on the noshear measurement
                if self.select is not None:
                    w=self._do_select(data['mcal_s2n'][wfield])
                    w=wfield[w]
                else:
                    w=wfield

                ntot  += wfield.size
                nkeep += w.size


                M1=data['mcal_pars'][w,2]
                M2=data['mcal_pars'][w,3]

                wt=self._get_noise_weight_mom(data, w, 'noshear')

                sums['wsum'][i] += wt.sum()
                sums['g'][i,0]  += (M1*wt).sum()
                sums['g'][i,1]  += (M2*wt).sum()

                if 'mcal_gpsf' in data.dtype.names:
                    sums['gpsf'][i] += data['mcal_gpsf'][w].sum(axis=0)

                # must use the same weight for these, equivalent to same
                # selection
                for type in ngmix.metacal.METACAL_TYPES:
                    if type=='noshear':
                        continue
                    mcalname='mcal_pars_%s' % type

                    if mcalname in data.dtype.names:
                        sumname='g_%s' % type

                        M1 = data[mcalname][w,2]
                        M2 = data[mcalname][w,3]
                        sums[sumname][i,0] += (M1*wt).sum()
                        sums[sumname][i,1] += (M2*wt).sum()
                    else:
                        #print("    skipping:",mcalname)
                        pass

                # now the selection terms, always needed since
                # we are weighting
                for type in ngmix.metacal.METACAL_TYPES:
                    if type=='noshear':
                        continue
                    ts2n_name='mcal_s2n_%s' % type

                    if ts2n_name in data.dtype.names:
                        wsumname = 's_wsum_%s' % type
                        sumname = 's_g_%s' % type

                        if self.select is not None:
                            w=self._do_select(data[ts2n_name][wfield])
                            w=wfield[w]
                        else:
                            w=wfield

                        M1=data['mcal_pars'][w,2]
                        M2=data['mcal_pars'][w,3]

                        wt=self._get_noise_weight_mom(data, w, type)

                        sums[wsumname][i]  += wt.sum()
                        sums[sumname][i,0] += (M1*wt).sum()
                        sums[sumname][i,1] += (M2*wt).sum()
                    else:
                        #print("    skipping:",ts2n_name)
                        pass

        if self.select is not None:
            self._print_frac(ntot,nkeep)
        return sums


    def _average_sums(self, sums):
        """
        divide by sum of weights and get g for each field

        Also average the responses over all data
        """

        # g averaged in each field
        g = sums['g'].copy()
        gpsf = sums['gpsf'].copy()

        winv = 1.0/sums['wsum']
        g[:,0]    *= winv
        g[:,1]    *= winv
        gpsf[:,0] *= winv
        gpsf[:,1] *= winv

        # responses averaged over all fields
        R = zeros(2)
        Rpsf = zeros(2)
        Rsel = zeros(2)
        Rsel_psf = zeros(2)

        factor = 1.0/(2.0*self.step)

        wsum=sums['wsum'].sum()

        g1p = sums['g_1p'][:,0].sum()/wsum
        g1m = sums['g_1m'][:,0].sum()/wsum
        g2p = sums['g_2p'][:,1].sum()/wsum
        g2m = sums['g_2m'][:,1].sum()/wsum

        g1p_psf = sums['g_1p_psf'][:,0].sum()/wsum
        g1m_psf = sums['g_1m_psf'][:,0].sum()/wsum
        g2p_psf = sums['g_2p_psf'][:,1].sum()/wsum
        g2m_psf = sums['g_2m_psf'][:,1].sum()/wsum

        R[0] = (g1p - g1m)*factor
        R[1] = (g2p - g2m)*factor
        Rpsf[0] = (g1p_psf - g1m_psf)*factor
        Rpsf[1] = (g2p_psf - g2m_psf)*factor

        #Rmean = R.mean()
        #R[:] = Rmean
        #R[0] = R[1]
        #R[:] = R[1]
        #R *= (-1)
        #R[:]=1

        print("R:",R)
        print("Rpsf:",Rpsf)

        # selection terms
        if self.do_selection:
            s_g1p = sums['s_g_1p'][:,0].sum()/sums['s_wsum_1p'].sum()
            s_g1m = sums['s_g_1m'][:,0].sum()/sums['s_wsum_1m'].sum()
            s_g2p = sums['s_g_2p'][:,1].sum()/sums['s_wsum_2p'].sum()
            s_g2m = sums['s_g_2m'][:,1].sum()/sums['s_wsum_2m'].sum()

            Rsel[0] = (s_g1p - s_g1m)*factor
            Rsel[1] = (s_g2p - s_g2m)*factor

            # can be zero if we aren't calculating psf terms (roundified psf)
            tsum=sums['s_wsum_1p_psf'].sum()
            if tsum != 0.0:
                s_g1p_psf = sums['s_g_1p_psf'][:,0].sum()/sums['s_wsum_1p_psf'].sum()
                s_g1m_psf = sums['s_g_1m_psf'][:,0].sum()/sums['s_wsum_1m_psf'].sum()
                s_g2p_psf = sums['s_g_2p_psf'][:,1].sum()/sums['s_wsum_2p_psf'].sum()
                s_g2m_psf = sums['s_g_2m_psf'][:,1].sum()/sums['s_wsum_2m_psf'].sum()

                Rsel_psf[0] = (s_g1p_psf - s_g1m_psf)*factor
                Rsel_psf[1] = (s_g2p_psf - s_g2m_psf)*factor

            print()
            print("Rsel:",Rsel)
            print("Rpsf_sel:",Rsel_psf)

        return g, gpsf, R, Rpsf, Rsel, Rsel_psf

    def _print_frac(self, ntot, nkeep):
        frac=float(nkeep)/ntot
        print("        kept: %d/%d = %g" % (nkeep,ntot,frac))

    def _do_select_old(self, s2n):
        """
        currently only s/n
        """
        logic=eval(self.select)
        w,=numpy.where(logic)
        return w

    def _do_select(self, data, w, type=None):
        """
        currently only s/n
        """

        s2n=self._get_s2n(data, w, type=type)
        size=self._get_size(data, w, type=type)

        if self.namer('flux_s2n') in data.dtype.names:
            flux_s2n = self._get_flux_s2n(data, w, type=type)

        logic=eval(self.select)
        w,=numpy.where(logic)
        return w

    def _get_s2n_name(self, data, type=None):
        n=self.namer
        if n('s2n_r') in data.dtype.names:
            name=n('s2n_r')
        elif n('s2n') in data.dtype.names:
            name=n('s2n')
        elif n('s2n_w') in data.dtype.names:
            print("Using s2n_w for selections")
            name=n('s2n_w')
        else:
            return None

        if type is not None and type != 'noshear':
            name='%s_%s' % (name, type)

        return name

    def _get_s2n(self, data, w, type=None):
        name=self._get_s2n_name(data, type=type)
        if name is None:
            return None

        return data[name][w]

    def _get_size(self, data, w, type=None):
        name=self._get_pars_name(data, type=type)
        return data[name][w,4]

    def _get_pars_name(self, data, type=None):
        n=self.namer
        name=n('pars')
        if type is not None and type != 'noshear':
            name='%s_%s' % (name, type)

        return name

    def _get_flux_s2n(self, data, w, type=None):
        name=self._get_flux_s2n_name(data, type=type)
        return data[name][w]


    def _get_flux_s2n_name(self, data, type=None):
        n=self.namer
        name=n('flux_s2n')
        if type is not None and type != 'noshear':
            name='%s_%s' % (name, type)

        return name


    def _read_means(self):
        fname=self._get_means_file()
        print("reading:",fname)
        return fitsio.read(fname)

    def _write_means(self):
        fname=self._get_means_file()
        eu.ostools.makedirs_fromfile(fname)
        #print("writing:",fname)
        fitsio.write(fname, self.means, clobber=True)

    def _get_fname_extra(self, last=None):
        runs=self.args.runs
        if len(runs) > 1:
            extra=runs[1:]
        else:
            extra=[]

        if self.args.weighted:
            extra += ['weighted']

        if self.args.select is not None:
            s=self.select.replace(' ','-').replace('(','').replace(')','').replace('[','').replace(']','').replace('"','').replace("'",'')
            extra += ['select',s]

        if self.args.ntest is not None:
            extra += ['test%d' % self.args.ntest]

        if last is not None:
            extra += [last]

        if len(extra) > 0:
            extra = '-'.join(extra)
        else:
            extra=None

        return extra



    def _get_means_file(self):

        extra=self._get_fname_extra()
        fname=files.get_means_url(self.args.runs[0], extra=extra)
        return fname


    def _get_fit_plot_file(self):
        extra=self._get_fname_extra(last='fit-m-c')
        fname=files.get_plot_url(self.args.runs[0], extra=extra)
        return fname

    def _get_resid_hist_file(self):
        extra=self._get_fname_extra(last='resid-hist')
        fname=files.get_plot_url(self.args.runs[0], extra=extra)
        return fname


    def _get_sums_struct(self):
        dt=self._get_sums_dt()
        return numpy.zeros(self['nshear'], dtype=dt)

    def _get_sums_dt(self):
        dt=[
            ('wsum','f8'),
            ('g','f8',2),
            ('gpsf','f8',2),

            ('g_1p','f8',2),
            ('g_1m','f8',2),
            ('g_2p','f8',2),
            ('g_2m','f8',2),
            ('g_1p_psf','f8',2),
            ('g_1m_psf','f8',2),
            ('g_2p_psf','f8',2),
            ('g_2m_psf','f8',2),

            # selection terms
            ('s_wsum_1p','f8'),
            ('s_wsum_1m','f8'),
            ('s_wsum_2p','f8'),
            ('s_wsum_2m','f8'),
            ('s_g_1p','f8',2),
            ('s_g_1m','f8',2),
            ('s_g_2p','f8',2),
            ('s_g_2m','f8',2),

            ('s_wsum_1p_psf','f8'),
            ('s_wsum_1m_psf','f8'),
            ('s_wsum_2p_psf','f8'),
            ('s_wsum_2m_psf','f8'),
            ('s_g_1p_psf','f8',2),
            ('s_g_1m_psf','f8',2),
            ('s_g_2p_psf','f8',2),
            ('s_g_2m_psf','f8',2),

            ('shear_true','f8',2),
        ]
        return dt

    def _set_select(self):
        self.select=None
        self.do_selection=False

        if self.args.select is not None:
            self.select = self.args.select
            self.do_selection=True
        elif self.args.select_from is not None:
            with open(self.args.select_from) as fobj:
                d=yaml.load(fobj)
            self.do_selection=True

            self.select = d['select'].strip()
        elif self.args.weighted:
            self.do_selection=True


    def plot_fits(self):
        import biggles
        biggles.configure('default','fontsize_min',1.5)

        means=self.means
        fits=self.fits
        args=self.args
        #Q=calc_q(fits)

        if args.yrange is not None:
            yrange=[float(r) for r in args.yrange.split(',')]
        else:
            yrange=[-0.01,0.01]

        xrng=args.xrange
        if xrng is not None:
            xrng=[float(r) for r in args.xrange.split(',')]

        tab=biggles.Table(1,2)
        tab.aspect_ratio=0.5

        diff = means['shear'] - means['shear_true']

        plts=[]
        for i in [0,1]:

            x = means['shear_true'][:,i]
            plt =biggles.plot(
                x,
                diff[:,i],
                xlabel=r'$\gamma_{%d}$ true' % (i+1,),
                ylabel=r'$\Delta \gamma_{%d}$' % (i+1,),
                yrange=yrange,
                xrange=xrng,
                visible=False,
            )
            yfit=fits['m'][0,i]*x + fits['c'][0,i]

            z=biggles.Curve(x, x*0, color='black')
            c=biggles.Curve(x, yfit, color='red')
            plt.add(z,c)

            '''
            mstr='m%d: %.2g +/- %.2g' % (i+1,fits['m'][0,i],fits['merr'][0,i])
            cstr='c%d: %.2g +/- %.2g' % (i+1,fits['c'][0,i],fits['cerr'][0,i])
            mlab=biggles.PlotLabel(0.1,0.9,
                                   mstr,
                                   halign='left')
            clab=biggles.PlotLabel(0.1,0.85,
                                   cstr,
                                   halign='left')
            plt.add(mlab,clab)
            '''
            if False and i==0:
                Qstr='Q: %d' % (int(Q),)
                Qlab=biggles.PlotLabel(0.1,0.8,
                                       Qstr,
                                       halign='left')
                plt.add(Qlab)


            tab[0,i] = plt

        fname=self._get_fit_plot_file()
        eu.ostools.makedirs_fromfile(fname)
        print("writing:",fname)
        tab.write_eps(fname)

        if args.show:
            tab.show(width=1000, height=1000)


    def plot_resid_hist(self):
        import biggles

        means=self.means
        fits=self.fits
        args=self.args
        #Q=calc_q(fits)

        diff = means['shear'] - means['shear_true']

        plt = biggles.plot_hist(diff, nbin=20, visible=False,
                               xlabel=r'$\gamma - \gamma_{True}$')

        dmax=numpy.abs(diff).max() 
        plt.xrange=[-1.3*dmax, 1.3*dmax]

        fname=self._get_resid_hist_file()
        eu.ostools.makedirs_fromfile(fname)
        print("writing:",fname)
        plt.write_eps(fname)

        if args.show:
            plt.show(width=1000, height=1000)

class SummerNSim(Summer):
    """
    For NSim we have the shears in the simulation config
    """
    def __init__(self, conf, args):

        shear_pdf = shearpdf.get_shear_pdf(conf['simc'])
        shears=shear_pdf.shears
        super(SummerNSim,self).__init__(conf, args, shears=shears)


class SummerMoments(SummerNSim):
    def go(self):

        if self.args.fit_only:
            self.means=self._read_means()
        else:

            extra=self._get_fname_extra()

            sums = self.do_sums()

            args=self.args

            e,R,Rsel=self._average_sums(sums)

            means=get_mean_struct(self['nshear'])
            means_nocorr=get_mean_struct(self['nshear'])

            for i in xrange(self['nshear']):

                shear_true = self.shears[i]

                emean = e[i]

                shear        = emean/(R+Rsel)
                shear_nocorr = emean/R

                means['shear'][i] = shear
                means['shear_err'][i] = 1.0
                if isinstance(shear_true,ngmix.Shape):
                    means['shear_true'][i,0] = shear_true.g1
                    means['shear_true'][i,1] = shear_true.g2
                else:
                    means['shear_true'][i] = shear_true

                means_nocorr['shear'][i] = shear_nocorr
                means_nocorr['shear_err'][i] = 1.0
                means_nocorr['shear_true'][i] = means['shear_true'][i]

            self.means=means
            self.means_nocorr=means_nocorr
            self._write_means()

        if self.do_selection:
            print("without correction")
            junk=reredux.averaging.fit_m_c(self.means_nocorr)
            junk=reredux.averaging.fit_m_c(self.means_nocorr,onem=True)
            print("\nwith correction")

        self.fits=reredux.averaging.fit_m_c(self.means)
        self.fitsone=reredux.averaging.fit_m_c(self.means,onem=True)


    def _get_M_T(self, data, w, type):
        if type=='noshear':
            name='mcal_pars'
        else:
            name='mcal_pars_%s' % type

        #Finv = 1.0/data[name][w,5]
        #Finv = 1.0/data['pars'][w,5]

        M = data[name][w,2:2+2].copy()
        T = data[name][w,4].copy()

        #M[:,0] *= Finv
        #M[:,1] *= Finv
        #T      *= Finv

        return M,T


    def _get_M_T_as_e(self, data, w, type):
        if type=='noshear':
            name='mcal_pars'
        else:
            name='mcal_pars_%s' % type

        #Finv = 1.0/data[name][w,5]
        #Finv = 1.0/data['pars'][w,5]

        M = data[name][w,2:2+2].copy()
        T = data[name][w,4].copy()

        M[:,0] /= T
        M[:,1] /= T

        T[:] = 1.0

        return M,T


    def do_sums1(self, data, sums=None):
        """
        just a binner and summer, no logic here
        """

        s2n_name='mcal_s2n'

        nshear=self['nshear']
        args=self.args

        h,rev = eu.stat.histogram(data['shear_index'],
                                  min=0,
                                  max=nshear-1,
                                  rev=True)
        nind = h.size
        assert nshear==nind

        if sums is None:
            sums=self._get_sums_struct()

        ntot=0
        nkeep=0
        for i in xrange(nshear):
            if rev[i] != rev[i+1]:
                wfield=rev[ rev[i]:rev[i+1] ]

                # first select on the noshear measurement
                if self.select is not None:
                    w=self._do_select(data[s2n_name][wfield])
                    w=wfield[w]
                else:
                    w=wfield

                ntot  += wfield.size
                nkeep += w.size

                sums['wsum'][i] += w.size

                M,T = self._get_M_T(data, w, 'noshear')
                sums['M'][i] += M.sum(axis=0)
                sums['T'][i] += T.sum()

                for type in ngmix.metacal.METACAL_TYPES:
                    if type=='noshear':
                        continue

                    Msumname='M_%s' % type
                    Tsumname='T_%s' % type

                    name='mcal_pars_%s' % type
                    if name in data.dtype.names:
                        M,T = self._get_M_T(data, w, type)
                        sums[Msumname][i] += M.sum(axis=0)
                        sums[Tsumname][i] += T.sum()

                # now the selection terms
                if self.select is not None:
                    for type in ngmix.metacal.METACAL_TYPES:
                        if type=='noshear':
                            continue

                        ts2n_name='%s_%s' % (s2n_name,type)

                        if ts2n_name in data.dtype.names:

                            wsumname = 's_wsum_%s' % type
                            Msumname = 's_M_%s' % type
                            Tsumname = 's_T_%s' % type

                            w=self._do_select(data[ts2n_name][wfield])
                            w=wfield[w]
                            sums[wsumname][i] += w.size

                            M,T = self._get_M_T(data, w, 'noshear')
                            sums[Msumname][i] += M.sum(axis=0)
                            sums[Tsumname][i] += T.sum()
                        else:
                            pass

        if self.select is not None:
            self._print_frac(ntot,nkeep)
        return sums




    def _average_sums(self, sums):
        """
        divide by sum of weights and get g for each field

        Also average the responses over all data
        """

        M = sums['M'].copy()
        T = sums['T'].copy()

        #
        # averaged in each field
        #


        winv = 1.0/sums['wsum']
        M[:,0] *= winv
        M[:,1] *= winv
        T      *= winv

        Ta = T[:,numpy.newaxis]
        e = 0.5*M/Ta

        # overall means for responses

        wsum=sums['wsum'].sum()
        M_mean = sums['M'].sum(axis=0)/wsum
        T_mean = sums['T'].sum()/wsum
        T_meana = zeros(2) + T_mean

        factor = 1.0/(2.0*self.step)

        M1p = sums['M_1p'][:,0].sum()/wsum
        M1m = sums['M_1m'][:,0].sum()/wsum
        M2p = sums['M_2p'][:,1].sum()/wsum
        M2m = sums['M_2m'][:,1].sum()/wsum

        T1p = sums['T_1p'].sum()/wsum
        T1m = sums['T_1m'].sum()/wsum
        T2p = sums['T_2p'].sum()/wsum
        T2m = sums['T_2m'].sum()/wsum

        RM = zeros(2)
        RT = zeros(2)

        RM[0] = (M1p - M1m)*factor
        RM[1] = (M2p - M2m)*factor

        RT[0] = (T1p - T1m)*factor
        RT[1] = (T2p - T2m)*factor

        R = (1.0/T_meana) * RM  -  (M_mean/T_meana**2) * RT
        R *= 0.5

        print("RM:",RM)
        print("RT:",RT)
        print("R: ",R)

        # selection terms
        if self.do_selection:

            RMsel = zeros(2)
            RTsel = zeros(2)

            s_M1p = sums['s_M_1p'][:,0].sum()/sums['s_wsum_1p'].sum()
            s_M1m = sums['s_M_1m'][:,0].sum()/sums['s_wsum_1m'].sum()
            s_M2p = sums['s_M_2p'][:,1].sum()/sums['s_wsum_2p'].sum()
            s_M2m = sums['s_M_2m'][:,1].sum()/sums['s_wsum_2m'].sum()

            s_T1p = sums['s_T_1p'].sum()/sums['s_wsum_1p'].sum()
            s_T1m = sums['s_T_1m'].sum()/sums['s_wsum_1m'].sum()
            s_T2p = sums['s_T_2p'].sum()/sums['s_wsum_2p'].sum()
            s_T2m = sums['s_T_2m'].sum()/sums['s_wsum_2m'].sum()


            RMsel[0] = (s_M1p - s_M1m)*factor
            RMsel[1] = (s_M2p - s_M2m)*factor
            RTsel[0] = (s_T1p - s_T1m)*factor
            RTsel[1] = (s_T2p - s_T2m)*factor

            Rsel = (1.0/T_meana) * RMsel  -  (M_mean/T_meana**2) * RTsel
            Rsel *= 0.5

            print("Rsel:",Rsel)
        else:
            Rsel=zeros(2)

        return e, R, Rsel

    def _get_sums_dt(self):
        dt=[
            ('wsum','f8'),
            ('M','f8',2),
            ('T','f8'),

            ('M_1p','f8',2),
            ('M_1m','f8',2),
            ('M_2p','f8',2),
            ('M_2m','f8',2),

            ('T_1p','f8'),
            ('T_1m','f8'),
            ('T_2p','f8'),
            ('T_2m','f8'),


            # selection terms
            ('s_wsum_1p','f8'),
            ('s_wsum_1m','f8'),
            ('s_wsum_2p','f8'),
            ('s_wsum_2m','f8'),
            ('s_M_1p','f8',2),
            ('s_M_1m','f8',2),
            ('s_M_2p','f8',2),
            ('s_M_2m','f8',2),

            ('s_T_1p','f8'),
            ('s_T_1m','f8'),
            ('s_T_2p','f8'),
            ('s_T_2m','f8'),

        ]
        return dt

class SummerMomentsNoNorm(SummerMoments):

    def _get_M(self, data, w, type):
        if type=='noshear':
            name='mcal_pars'
        else:
            name='mcal_pars_%s' % type

        M = data[name][w,2:2+2].copy()
        return M

    def do_sums1(self, data, sums=None):
        """
        just a binner and summer, no logic here
        """

        s2n_name='mcal_s2n'

        nshear=self['nshear']
        args=self.args

        h,rev = eu.stat.histogram(data['shear_index'],
                                  min=0,
                                  max=nshear-1,
                                  rev=True)
        nind = h.size
        assert nshear==nind

        if sums is None:
            sums=self._get_sums_struct()

        ntot=0
        nkeep=0
        for i in xrange(nshear):
            if rev[i] != rev[i+1]:
                wfield=rev[ rev[i]:rev[i+1] ]

                # first select on the noshear measurement
                if self.select is not None:
                    w=self._do_select(data[s2n_name][wfield])
                    w=wfield[w]
                else:
                    w=wfield

                ntot  += wfield.size
                nkeep += w.size

                sums['wsum'][i] += w.size

                M = self._get_M(data, w, 'noshear')
                sums['M'][i] += M.sum(axis=0)

                for type in ngmix.metacal.METACAL_TYPES:
                    if type=='noshear':
                        continue

                    Msumname='M_%s' % type

                    name='mcal_pars_%s' % type
                    if name in data.dtype.names:
                        M = self._get_M(data, w, type)
                        sums[Msumname][i] += M.sum(axis=0)

                # now the selection terms
                if self.select is not None:
                    for type in ngmix.metacal.METACAL_TYPES:
                        if type=='noshear':
                            continue

                        ts2n_name='%s_%s' % (s2n_name,type)

                        if ts2n_name in data.dtype.names:

                            wsumname = 's_wsum_%s' % type
                            Msumname = 's_M_%s' % type

                            w=self._do_select(data[ts2n_name][wfield])
                            w=wfield[w]
                            sums[wsumname][i] += w.size

                            M = self._get_M(data, w, 'noshear')
                            sums[Msumname][i] += M.sum(axis=0)
                        else:
                            pass

        if self.select is not None:
            self._print_frac(ntot,nkeep)
        return sums




    def _average_sums(self, sums):
        """
        divide by sum of weights and get g for each field

        Also average the responses over all data
        """

        M = sums['M'].copy()

        #
        # averaged in each field
        #


        winv = 1.0/sums['wsum']
        M[:,0] *= winv
        M[:,1] *= winv

        wsum=sums['wsum'].sum()

        factor = 1.0/(2.0*self.step)

        M1p = sums['M_1p'][:,0].sum()/wsum
        M1m = sums['M_1m'][:,0].sum()/wsum
        M2p = sums['M_2p'][:,1].sum()/wsum
        M2m = sums['M_2m'][:,1].sum()/wsum

        R = zeros(2)

        R[0] = (M1p - M1m)*factor
        R[1] = (M2p - M2m)*factor

        print("R: ",R)

        # selection terms
        if self.do_selection:

            Rsel = zeros(2)

            s_M1p = sums['s_M_1p'][:,0].sum()/sums['s_wsum_1p'].sum()
            s_M1m = sums['s_M_1m'][:,0].sum()/sums['s_wsum_1m'].sum()
            s_M2p = sums['s_M_2p'][:,1].sum()/sums['s_wsum_2p'].sum()
            s_M2m = sums['s_M_2m'][:,1].sum()/sums['s_wsum_2m'].sum()

            Rsel[0] = (s_M1p - s_M1m)*factor
            Rsel[1] = (s_M2p - s_M2m)*factor

            print("Rsel:",Rsel)
        else:
            Rsel=zeros(2)

        return M, R, Rsel

    def _get_sums_dt(self):
        dt=[
            ('wsum','f8'),
            ('M','f8',2),

            ('M_1p','f8',2),
            ('M_1m','f8',2),
            ('M_2p','f8',2),
            ('M_2m','f8',2),

            # selection terms
            ('s_wsum_1p','f8'),
            ('s_wsum_1m','f8'),
            ('s_wsum_2p','f8'),
            ('s_wsum_2m','f8'),
            ('s_M_1p','f8',2),
            ('s_M_1m','f8',2),
            ('s_M_2p','f8',2),
            ('s_M_2m','f8',2),

        ]
        return dt


class AMSummer(SummerNSim):
    doM=False

    def _set_select(self):
        super(AMSummer,self)._set_select()
        #if self.doM:
        #    return

        select=self.select
        if select is None:
            select=[]
        else:
            select = ['( ' + select +' )']

        select += ['(T > 0.01)']

        self.select = ' & '.join(select)
        print("selection:",self.select)

        self.do_selection=True

    def _do_select(self, data, w, type):
        """
        currently only s/n
        """
        from numpy import abs

        s2n_name='mcal_s2n'
        T_name = 'mcal_T'

        if type != 'noshear':
            s2n_name = '%s_%s' % (s2n_name, type)
            T_name   = '%s_%s' % (T_name,   type)

        s2n = data[s2n_name][w]
        T   = data[T_name][w]

        g=self._get_g(data, w, type)
        e1 = g[:,0]*2
        e2 = g[:,1]*2

        e = numpy.sqrt(e1**2 + e2**2)

        logic=eval(self.select)

        wsel,=numpy.where(logic)

        return w[wsel]

    def _get_gpsf(self, data, w):
        Tpsf = data['mcal_psf_icc'][w]+data['mcal_psf_irr'][w]
        Tpsf_inv=1.0/Tpsf

        e1psf = (data['mcal_psf_icc'][w]-data['mcal_psf_irr'][w])*Tpsf_inv
        e2psf = 2*data['mcal_psf_irc'][w]*Tpsf_inv

        g1psf, g2psf = ngmix.shape.e1e2_to_g1g2(e1psf, e2psf)

        gpsf=numpy.zeros( (w.size, 2))

        gpsf[:,0] = g1psf
        gpsf[:,1] = g2psf

        return gpsf

    def _get_gpsf_frome(self, data, w):
        Tpsf = data['mcal_psf_icc'][w]+data['mcal_psf_irr'][w]
        Tpsf_inv=1.0/Tpsf

        g1psf = 0.5*(data['mcal_psf_icc'][w]-data['mcal_psf_irr'][w])*Tpsf_inv
        g2psf = 0.5*2*data['mcal_psf_irc'][w]*Tpsf_inv

        gpsf=numpy.zeros( (w.size, 2))

        gpsf[:,0] = g1psf
        gpsf[:,1] = g2psf

        return gpsf

    def _get_M(self, data, w, type):
        irr_name='mcal_irr'
        irc_name='mcal_irc'
        icc_name='mcal_icc'
        if type != 'noshear':
            irr_name='%s_%s' % (irr_name,type)
            irc_name='%s_%s' % (irc_name,type)
            icc_name='%s_%s' % (icc_name,type)

        irr = data[irr_name][w]
        irc = data[irc_name][w]
        icc = data[icc_name][w]

        M1 = icc - irr
        M2 = 2.0*irc
        T= irr + icc

        return M1, M2, T


    def _get_e(self, data, w, type):
        """
        really getting e, but multiply by 0.5 to
        approximately get in right scale as g
        """

        M1, M2, T = self._get_M(data, w, type)

        Tinv=1.0/T

        # 0.5 to approximately put in g units, just so we can
        # more easily examine the response
        e=numpy.zeros( (w.size, 2) )
        e[:,0] = M1*Tinv
        e[:,1] = M2*Tinv

        return e

    def _get_g(self, data, w, type):
        g=self._get_e(data, w,type)
        if g is not None:
            g *= 0.5
        return g

    def do_sums1(self, data, sums=None):
        """
        just a binner and summer, no logic here
        """


        s2n_name='mcal_s2n'
        T_name = 'mcal_T'

        nshear=self['nshear']
        args=self.args

        h,rev = eu.stat.histogram(data['shear_index'],
                                  min=0,
                                  max=nshear-1,
                                  rev=True)
        nind = h.size
        assert nshear==nind

        if sums is None:
            sums=self._get_sums_struct()

        ntot=0
        nkeep=0
        for i in xrange(nshear):
            if rev[i] != rev[i+1]:
                wfield=rev[ rev[i]:rev[i+1] ]

                # first select on the noshear measurement
                if self.select is not None:
                    w=self._do_select(data, wfield, 'noshear')
                else:
                    w=wfield

                ntot  += wfield.size
                nkeep += w.size

                sums['wsum'][i] += w.size

                if self.doM:
                    M1,M2,T = self._get_M(data, w, 'noshear')
                    sums['g'][i,0]    += M1.sum()
                    sums['g'][i,1]    += M2.sum()
                else:
                    g    = self._get_g(data, w, 'noshear')
                    sums['g'][i]    += g.sum(axis=0)

                gpsf = self._get_gpsf(data, w)
                sums['gpsf'][i] += gpsf.sum(axis=0)

                for type in ngmix.metacal.METACAL_TYPES:
                    if type=='noshear':
                        continue

                    sumname='g_%s' % type

                    if self.doM:
                        M1,M2,T = self._get_M(data, w, type)
                        sums[sumname][i,0] += M1.sum()
                        sums[sumname][i,1] += M2.sum()
                    else:
                        g = self._get_g(data, w, type)
                        if g is not None:
                            sums[sumname][i] += g.sum(axis=0)

                # now the selection terms
                if self.select is not None:
                    for type in ngmix.metacal.METACAL_TYPES:
                        if type=='noshear':
                            continue

                        ts2n_name='%s_%s' % (s2n_name,type)

                        if ts2n_name in data.dtype.names:

                            wsumname = 's_wsum_%s' % type
                            sumname = 's_g_%s' % type

                            w=self._do_select(data, wfield, type)
                            sums[wsumname][i] += w.size

                            if self.doM:
                                M1,M2,T = self._get_M(data, w, 'noshear')
                                sums[sumname][i,0] += M1.sum()
                                sums[sumname][i,1] += M2.sum()
                            else:
                                g = self._get_g(data, w, 'noshear')
                                sums[sumname][i] += g.sum(axis=0)
                        else:
                            pass

        if self.select is not None:
            self._print_frac(ntot,nkeep)
        return sums


# quick line fit pulled from great3-public code
def _calculateSvalues(xarr, yarr, sigma2=1.):
    """Calculates the intermediate S values required for basic linear regression.

    See, e.g., Numerical Recipes (Press et al 1992) Section 15.2.
    """
    if len(xarr) != len(yarr):
        raise ValueError("Input xarr and yarr differ in length!")
    if len(xarr) <= 1:
        raise ValueError("Input arrays must have 2 or more values elements.")

    S = len(xarr) / sigma2
    Sx = numpy.sum(xarr / sigma2)
    Sy = numpy.sum(yarr / sigma2)
    Sxx = numpy.sum(xarr * xarr / sigma2)
    Sxy = numpy.sum(xarr * yarr / sigma2)
    return (S, Sx, Sy, Sxx, Sxy)

def fitline(xarr, yarr):
    """Fit a line y = a + b * x to input x and y arrays by least squares.

    Returns the tuple (a, b, Var(a), Cov(a, b), Var(b)), after performing an internal estimate of
    measurement errors from the best-fitting model residuals.

    See Numerical Recipes (Press et al 1992; Section 15.2) for a clear description of the details
    of this simple regression.
    """
    # Get the S values (use default sigma2, best fit a and b still valid for stationary data)
    S, Sx, Sy, Sxx, Sxy = _calculateSvalues(xarr, yarr)
    # Get the best fit a and b
    Del = S * Sxx - Sx * Sx
    a = (Sxx * Sy - Sx * Sxy) / Del
    b = (S * Sxy - Sx * Sy) / Del
    # Use these to estimate the sigma^2 by residuals from the best-fitting model
    ymodel = a + b * xarr
    sigma2 = numpy.mean((yarr - ymodel)**2)
    # And use this to get model parameter error estimates
    var_a  = sigma2 * Sxx / Del
    cov_ab = - sigma2 * Sx / Del
    var_b  = sigma2 * S / Del

    a_err = numpy.sqrt(var_a)
    b_err = numpy.sqrt(var_b)

    return {'offset':a,
            'offset_err':a_err,
            'slope':b,
            'slope_err':b_err,
            'cov':cov_ab}
    #return a, a_err, b, b_err, cov_ab

def fitline_zero_offset(x, y):

    # Our model is y = a * x, so things are quite simple, in this case...
    # x needs to be a column vector instead of a 1D vector for this, however.
    x = x[:,numpy.newaxis]
    a, _, _, _ = numpy.linalg.lstsq(x, y)

    return {'offset':0.0,
            'offset_err':0.0,
            'slope':a[0],
            'slope_err':0.0}

def plot_line_fit(args, extra, x, y, res, xlabel, ylabel, label_error=True):
    import biggles
    plt=biggles.FramedPlot()

    ymin=y.min()
    ymax=y.max()
    if ymin < 0:
        yr = [1.1*ymin, 0.0]
    else:
        yr = [0, 1.1*ymax]

    xr = [0.0, 1.1*x.max()]

    plt.xrange=xr
    plt.yrange=yr
    plt.xlabel=xlabel
    plt.ylabel=ylabel
    plt.aspect_ratio=1

    xfit = numpy.linspace(0, xr[1])
    yfit = res['offset'] + res['slope']*xfit

    pts = biggles.Points(x,y,type='filled circle')
    c = biggles.Curve(xfit, yfit, color='blue')

    if label_error:
        alab=r'$slope = %.3g \pm %.3g' % (res['slope'],res['slope_err'])
        blab=r'$offset = %.3g \pm %.3g' % (res['offset'],res['offset_err'])
    else:
        alab=r'$slope = %.3g' % (res['slope'],)
        blab=r'$offset = %.3g' % (res['offset'],)
    alabel=biggles.PlotLabel(0.9, 0.9, alab, halign='right')
    blabel=biggles.PlotLabel(0.9, 0.85, blab, halign='right')

    plt.add(c, pts, alabel, blabel)

    plotfile=files.get_plot_url(args.runs[0], extra)

    print("writing:",plotfile)
    eu.ostools.makedirs_fromfile(plotfile)
    plt.write_eps(plotfile)

    if args.show:
        plt.show()

def print_Rs(R, Rerr, Rpsf, Rpsf_err, type=''):
    p='%s: %.5g +/- %.5g'
    for i in xrange(2):
        n='R%s_psf[%d]' % (type,i+1)
        print(p % (n,Rpsf[i],Rpsf_err[i]))

    for i in xrange(2):
        for j in xrange(2):
            n='R%s[%d,%d]' % (type,(i+1),(j+1))
            print(p % (n,R[i,j],Rerr[i,j]))

def get_s2n_weights_old(s2n, args):
    #return numpy.ones(s2n.size)
    #print("s2n soft:",args.s2n_soft)
    wts = 1.0/(1.0 + (args.s2n_soft/s2n)**2 )
    return wts

def get_s2n_weights(s2n, args):
    #return numpy.ones(s2n.size)
    #print("s2n soft:",args.s2n_soft)
    sigma=5.0
    x=zeros(s2n.size)
    wts=zeros(s2n.size)
    x[:] = s2n
    x -= 10.0
    x *= (1.0/sigma)
    ngmix._gmix.erf_array(x, wts)

    wts += 1.0
    wts *= 0.5
    return wts


def get_noise_weights(g_cov, args):
    wts = 1.0/(2*args.shapenoise**2 + g_cov[:,0,0] + g_cov[:,1,1])

    '''
    w,=numpy.where( numpy.isnan(wts) )
    if w.size > 0:
        print("fixing %d/%d isnan" % (w.size, g_cov.shape[0]))
        wts[w] = 0.0
    '''
    return wts


def get_mean_struct(n):
    dt=[('shear','f8',2),
        ('shear_true','f8',2),
        ('shear_err','f8',2)]

    means = numpy.zeros(n, dtype=dt)
    return means

def get_boot_struct(nboot):
    dt=[('m','f8',2),
        ('c','f8',2)]

    bs = numpy.zeros(nboot, dtype=dt)
    return bs


