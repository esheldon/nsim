from __future__ import print_function
import sys
import os
from pprint import pprint
import numpy
from numpy import sqrt, array, diag, ones, zeros
from numpy import where, newaxis, exp, log
from numpy import newaxis
import yaml
import biggles

import fitsio

import ngmix
import nsim
from . import files

import argparse
import esutil as eu
from esutil.numpy_util import between

CHUNKSIZE=1000000

try:
    import reredux
except ImportError:
    pass

class Summer(dict):
    def __init__(self, conf, shears, args):
        self.update(conf)
        self.args=args

        self._set_select()

        self.step = self['metacal_pars'].get('step',0.01)

        self.shears = shears
        self['nshear']=len(self.shears)

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

            for i in xrange(self['nshear']):

                shear_true = self.shears[i]

                gmean = g[i]

                c        = (Rpsf+Rsel_psf)*gpsf[i]
                c_nocorr = Rpsf*gpsf[i]

                shear        = (gmean-c)/(R+Rsel)
                shear_nocorr = (gmean-c_nocorr)/R

                means['shear'][i] = shear
                means['shear_err'][i] = 1.0
                means['shear_true'][i] = shear_true

                means_nocorr['shear'][i] = shear_nocorr
                means_nocorr['shear_err'][i] = 1.0
                means_nocorr['shear_true'][i] = shear_true

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
        return fname

    def do_sums(self):

        args=self.args

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
                nchunks = nrows//CHUNKSIZE

                if (nrows % CHUNKSIZE) > 0:
                    nchunks += 1

                beg=0
                for i in xrange(nchunks):
                    print("    chunk %d/%d" % (i+1,nchunks))

                    end=beg+CHUNKSIZE

                    data = hdu[beg:end]

                    data=self._preselect(data)

                    ntot += data.size

                    #sums=self.do_sums1(data, sums=sums)
                    if True and 'mcal_g' not in data.dtype.names:
                        if False:
                            sums=self.do_sums1_moms_wt(data, sums=sums)
                        elif True:
                            sums=self.do_sums1_moms(data, sums=sums)
                        elif False:
                            sums=self.do_sums1_moms_psfcorr(data, sums=sums)
                    else:
                        sums=self.do_sums1(data, sums=sums)

                    beg = beg + CHUNKSIZE

                    if args.ntest is not None and ntot > args.ntest:
                        break

        return sums

    def _preselect(self, data):
        """
        sub-classes might make a pre-selection, e.g. of some flags
        """
        return data

        norig=data.size

        '''
        w,=where(
            (data['mcal_pars_cov'][:,0,0] > 0) &
            (data['mcal_pars_cov'][:,1,1] > 0) &
            (data['mcal_pars_cov'][:,2,2] > 0) &
            (data['mcal_pars_cov'][:,3,3] > 0) &
            (data['mcal_pars_cov'][:,4,4] > 0) &
            (data['mcal_pars_cov'][:,5,5] > 0)
        )
        '''

        #movsig1=data['mcal_pars'][:,0]/sqrt(data['mcal_pars_cov'][:,0,0])
        #movsig2=data['mcal_pars'][:,1]/sqrt(data['mcal_pars_cov'][:,1,1])

        w0,=where(data['mcal_pars'][:,4] > 0)

        e1=data['mcal_pars'][w0,2]/data['mcal_pars'][w0,4]
        e2=data['mcal_pars'][w0,3]/data['mcal_pars'][w0,4]

        w,=where(  (numpy.abs(e1) < 0.9999)
                 & (numpy.abs(e2) < 0.9999) )

        w=w0[w]
        w,=where(
            #(numpy.abs(movsig1) < 2 ) &
            #(numpy.abs(movsig2) < 2 ) &
            (numpy.abs(data['mcal_pars'][:,0]) < 2.5)  # would change with s/n
            (numpy.abs(data['mcal_pars'][:,1]) < 2.5 ) &
            (numpy.abs(data['mcal_pars'][:,2]) < 20.0) &
            (numpy.abs(data['mcal_pars'][:,3]) < 20.0) &
            (data['mcal_pars'][:,4] < 150.0)
        )

        print("kept %d/%d preselect" % (w.size, norig))

        data=data[w]
        return data

    def _get_s2n_name(self, data):
        if 'mcal_s2n_r' in data.dtype.names:
            s2n_name='mcal_s2n_r'
        elif 'mcal_s2n' in data.dtype.names:
            s2n_name='mcal_s2n'

    def _get_bname_and_beg(self, data):
        if 'mcal_g' in data.dtype.names:
            bname='mcal_g'
            beg=0
        else:
            bname='mcal_pars'
            beg=2

        return bname, beg

    def _get_g(self, data, w, type):
        bname, beg=self._get_bname_and_beg(data):
        if type=='noshear':
            name=bname
        else:
            name='%s_%s' % (bname, type)

        if name not in data.dtype.names:
            g=None
        else:
            g = data[name][w,beg:beg+2]
        return g

    def do_sums1(self, data, sums=None):
        """
        just a binner and summer, no logic here
        """

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
                g=self._get_g(data, w, 'noshear')
                sums['g'][i] += g.sum(axis=0)

                if 'mcal_gpsf' in data.dtype.names:
                    sums['gpsf'][i] += data['mcal_gpsf'][w].sum(axis=0)

                for type in ngmix.metacal.METACAL_TYPES:
                    if type=='noshear':
                        continue
                    mcalname='%s_%s' % (bname,type)

                    g=self._get_g(data, w, type)
                    if g is not None:
                        sumname='g_%s' % type

                        sums[sumname][i] += g.sum(axis=0)
                    else:
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

                            g=self._get_g(data, w, 'noshear')
                            sums[sumname][i] += g.sum(axis=0)
                        else:
                            #print("    skipping:",ts2n_name)
                            pass

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

    def _do_select(self, s2n):
        """
        currently only s/n
        """
        logic=eval(self.select)
        w,=numpy.where(logic)
        #print("   kept: %d/%d" % (w.size,s2n.size))
        return w

    def _read_means(self):
        fname=self._get_means_file()
        print("reading:",fname)
        return fitsio.read(fname)

    def _write_means(self):
        fname=self._get_means_file()
        eu.ostools.makedirs_fromfile(fname)
        print("writing:",fname)
        fitsio.write(fname, self.means, clobber=True)

    def _get_fname_extra(self, last=None):
        runs=self.args.runs
        if len(runs) > 1:
            extra=runs[1:]
        else:
            extra=[]

        if self.args.weights is not None:
            extra += ['wts', self.args.weights]

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

            # these get filled in at the end
            ('R','f8',2),
            ('Rpsf','f8',2),
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
        shears = conf['simc']['shear']['shears']

        super(SummerNSim,self).__init__(conf, shears, args)


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
    #return numpy.ones(s2n.size)
    #print("s2n soft:",args.s2n_soft)
    wts = 1.0/(2*args.shapenoise**2 + g_cov[:,0,0] + g_cov[:,1,1])
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


