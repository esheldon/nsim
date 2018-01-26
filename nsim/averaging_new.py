from __future__ import print_function
try:
    xrange
except:
    xrange=range

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
from . import files
from . import shearpdf

from . import util
from .util import Namer

import argparse
import esutil as eu
from esutil.numpy_util import between

class Summer(dict):
    def __init__(self, args):

        self.args=args
        self._load_config()
        self._load_shear_pdf()

        self.chunksize=args.chunksize

        self._set_select()

        self.gpsf_name='mcal_psfrec_g'
        self.gpsf_orig_name='psfrec_g'

        if args.R is not None:
            self.Rinput = array([float(R) for R in args.R.split(',')])
        if args.Rselect is not None:
            self.Rselect_input = array([float(R) for R in args.Rselect.split(',')])

    def _load_config(self):
        args=self.args
        if 'runs' in args.runs:
            # this is a runs config file
            data=files.read_config(args.runs)
            self.runs=data['runs']

        else:
            self.runs=[args.runs]

        conf = files.read_config(self.runs[0])
        conf['simc'] = files.read_config(conf['sim'])

        conf['simc']['do_ring'] = conf['simc'].get('do_ring',False)
        print("do_ring:",conf['simc']['do_ring'])

        self.update(conf)
        self._set_step()

    def _set_step(self):
        self.step = self['metacal_pars'].get('step',0.01)
        self.stepfac=1.0/(2*self.step)

    def _load_shear_pdf(self):
        if 'shear' in self['simc']:
            shear_pdf = shearpdf.get_shear_pdf(
                self['simc']['shear'],
                None,
            )
            self.shears=shear_pdf.shears
            self['nshear']=len(self.shears)
        else:
            self.shears=None
            # one shear, since we will take the average
            self['nshear'] = 1

    def get_namer(self, type):
        return Namer(front='mcal', back=type)

    def get_name(self, name, type):
        n=self.get_namer(type)
        return n(name)

    def go(self):

        if self.args.fit_only:
            self.means=self._read_means()
        else:

            sums = self.do_sums()

            args=self.args

            g,gerr,gpsf,gpsf_orig,R,Rpsf,Rsel,Rsel_psf=self._average_sums(sums)

            means=get_mean_struct(self['nshear'])
            means_nocorr=get_mean_struct(self['nshear'])
            
            wkeep,=numpy.where(sums['wsum'] > 0)
            print("total wsum:",sums['wsum'].sum())

            for i in xrange(self['nshear']):

                if self.shears is not None:
                    shear_true = self.shears[i]
                else:
                    shear_true=zeros(2)
                    shear_true[0] = sums['shear_sum'][i,0]/sums['R11_sum'][i]
                    shear_true[1] = sums['shear_sum'][i,1]/sums['R22_sum'][i]
                    print("shear_true:",shear_true)

                gmean = g[i]
                gmean_err = gerr[i]

                c        = (Rpsf+Rsel_psf)*gpsf[i]

                if self.args.corr_psf_orig:
                    print("correcting for orig psf",gpsf_orig[i])
                    c += Rsel_psf*gpsf_orig[i]

                c_nocorr = Rpsf*gpsf[i]

                shear        = (gmean-c)/(R+Rsel)
                shear_nocorr = (gmean-c_nocorr)/R


                shear_err = gmean_err/(R+Rsel)
                shear_err_nocorr = gmean_err/R

                print("shear_meas:",shear)
                means['shear'][i] = shear
                means['shear_err'][i] = shear_err
                if isinstance(shear_true,ngmix.Shape):
                    means['shear_true'][i,0] = shear_true.g1
                    means['shear_true'][i,1] = shear_true.g2
                else:
                    means['shear_true'][i] = shear_true

                means_nocorr['shear'][i] = shear_nocorr
                means_nocorr['shear_err'][i] = shear_err_nocorr
                means_nocorr['shear_true'][i] = means['shear_true'][i]

            means=means[wkeep]
            means_nocorr=means_nocorr[wkeep]
            self.means=means
            self.means_nocorr=means_nocorr
            self._write_means()

        if means.size == 1:
            if self.do_selection:
                print("without correction")
                junk=get_m_c_oneshear(self.means_nocorr,nsig=args.nsigma)
                print("\nwith correction")

            self.fits=get_m_c_oneshear(self.means,nsig=args.nsigma)

        else:
            if self.do_selection:
                print("without correction")
                junk=fit_m_c(self.means_nocorr)
                junk=fit_m_c(self.means_nocorr,onem=True)
                print("\nwith correction")

            self.fits=fit_m_c(self.means,nsig=args.nsigma)
            self.fitsone=fit_m_c(self.means,onem=True,nsig=args.nsigma)


    def get_run_output(self, run):
        """
        collated file
        """

        fname = files.get_output_url(run, itrial=self.args.index)
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

        sums=None
        ntot=0
        for run in self.runs:
            if args.ntest is not None and ntot > args.ntest:
                break

            run_sums=self._try_read_sums(run)

            if args.force or run_sums is None:
                # no cache found, we need to do the sums

                run_sums=self._get_sums_struct()

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

                        if args.index is not None:
                            data = self._add_true_shear(data)

                        if 'shear_index' not in data.dtype.names:
                            data=self._add_shear_index(data)
                        else:
                            w,=numpy.where(data['shear_index']==-1)
                            if w.size > 0:
                                data['shear_index'][w]=0

                        data=self._preselect(data)

                        ntot += data.size

                        run_sums=self.do_sums1(data, sums=run_sums)

                        beg = beg + chunksize

                        if args.ntest is not None and ntot > args.ntest:
                            break

                self._write_sums(run, run_sums)

            if sums is None:
                sums=run_sums.copy()
            else:
                sums=add_sums(sums, run_sums)

        return sums

    def _add_shear_index(self, data):
        if self['nshear'] != 1:
            raise ValueError("need nshear==1 to fake "
                             "the shear index")

        add_dt=[('shear_index','i2')]
        new_data = eu.numpy_util.add_fields(
            data,
            add_dt,
        )
        new_data['shear_index'] = 0
        return new_data

    def _try_read_sums(self, run):
        sums_file=self._get_sums_file(run)

        if not os.path.exists(sums_file):
            sums_file=self._get_sums_file(run, use_select_string=True)
            if not os.path.exists(sums_file):
                print("sums file not found:",sums_file)
                return None

        print("reading cached sums file:",sums_file)
        sums = fitsio.read(sums_file)
        return sums

    def _preselect(self, data):
        """
        sub-classes might make a pre-selection, e.g. of some flags
        """
        
        if self.args.preselect:
            w,=numpy.where(data['s2n_true'] > 5)
            print("kept %d/%d in preselect" % (w.size, data.size))
            data=data[w]
 
        return data


    def _get_weights(self, data, w, type):

        if self.args.weighted:

            g_cov_name=self.get_name('g_cov',type)
            g_cov=data[g_cov_name][w]
            wts=get_noise_weights(g_cov, self.args)

        else:
            wts=numpy.ones(w.size)

        wa = wts[:,numpy.newaxis]
        return wts, wa

    def _get_g(self, data, w, type):
        name = self.get_name('g', type)

        if name not in data.dtype.names:
            g = None
        else:
            g = data[name][w]

        return g

    def _get_gpsf(self, data, w):
        name=self.gpsf_name
        if name in data.dtype.names:
            gpsf=data[name][w]
        else:
            gpsf=None

        name=self.gpsf_orig_name
        if name in data.dtype.names:
            gpsf_orig=data[name][w]
        else:
            gpsf_orig=None

        return gpsf, gpsf_orig

    def do_sums1(self, data, sums=None):
        """
        just a binner and summer, no logic here
        """

        names=data.dtype.names

        if 'shear_true' in names or 'shear' in names:
            sumshear=True
            if 'shear_true' in names:
                nm = 'shear_true'
            elif 'shear' in names:
                nm = 'shear'
        else:
            sumshear=False


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

                if sumshear:
                    R11 = (data['mcal_g_1p'][w,0] - 
                           data['mcal_g_1m'][w,0])*self.stepfac
                    R22 = (data['mcal_g_2p'][w,1] - 
                           data['mcal_g_2m'][w,1])*self.stepfac

                    sums['shear_sum'][i,0] += (R11*data[nm][w,0]).sum()
                    sums['shear_sum'][i,1] += (R22*data[nm][w,1]).sum()
                    sums['R11_sum'][i] += R11.sum()
                    sums['R22_sum'][i] += R22.sum()

                g = self._get_g(data, w, 'noshear')
                wts, wa = self._get_weights(data, w, 'noshear')

                wttot += wts.sum()
                twtmax = wts.max()
                if twtmax > wtmax:
                    wtmax = twtmax

                sums['g'][i]    += (g*wa).sum(axis=0)
                sums['gsq'][i]  += (g**2*wa**2).sum(axis=0)
                sums['wsq'][i]  += (wa**2).sum(axis=0)

                sums['wsum'][i] += wts.sum()

                gpsf, gpsf_orig =self._get_gpsf(data, w)
                if gpsf is not None:
                    sums['gpsf'][i] += (gpsf*wa).sum(axis=0)
                if gpsf_orig is not None:
                    sums['gpsf_orig'][i] += (gpsf_orig*wa).sum(axis=0)

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
                if self.do_selection:

                    for type in ngmix.metacal.METACAL_TYPES:
                        if type=='noshear':
                            continue

                        g_name=self.get_name('g', type)

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





    def _average_sums(self, sums):
        """
        divide by sum of weights and get g for each field

        Also average the responses over all data
        """

        # g averaged in each field
        g = sums['g'].copy()
        gpsf = sums['gpsf'].copy()

        gsq = sums['gsq'].copy()
        wsq = sums['wsq'].copy()


        winv = 1.0/sums['wsum']
        wainv = winv[:,newaxis]

        g[:,0]    *= winv
        g[:,1]    *= winv
        gpsf[:,0] *= winv
        gpsf[:,1] *= winv

        if 'gpsf_orig' in sums.dtype.names:
            gpsf_orig = sums['gpsf_orig'].copy()
            gpsf_orig[:,0] *= winv
            gpsf_orig[:,1] *= winv
        else:
            gpsf_orig=None

        # sum(w*2g*2
        gerrsq_sum = gsq - g**2*wsq
        gerr = sqrt(gerrsq_sum)*wainv

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

        if self.args.etype:
            Rpsf *= 0.5
            #Rpsf *= (0.01/0.019998000199980003)
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

                if self.args.etype:
                    Rsel_psf *= 0.5
                    #Rsel_psf *= (0.01/0.019998000199980003)

            print()
            print("Rsel:",Rsel)
            print("Rpsf_sel:",Rsel_psf)

        if self.args.R is not None:
            print("ignoring calculated R, using input")
            R = self.Rinput
            Rsel = self.Rselect_input
            print("R:",R)
            print("Rsel:",Rsel)

        return g, gerr, gpsf, gpsf_orig, R, Rpsf, Rsel, Rsel_psf

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

        n=self.get_namer(type)

        s2n=self._get_s2n(n, data, w)
        size=self._get_size(n, data, w)
        flux=self._get_flux(n, data, w)
        am_flux=self._get_am_flux(n, data, w)
        am_flux_s2n=self._get_am_flux_s2n(n, data, w)
        flux_true=self._get_flux_true(data, w)

        T = self._get_T(n, data, w)

        Tpsf = self._get_psf_T(data, w)

        if Tpsf is not None:
            Tratio = T/Tpsf

        if n('flux_s2n') in data.dtype.names:
            flux_s2n = self._get_flux_s2n(n, data, w)

        logic=eval(self.select)
        if self['simc']['do_ring']:
            w=util.ring_select(logic)
        else:
            w,=numpy.where(logic)
        return w

    def _get_s2n_name(self, n, data):

        if n('s2n_r') in data.dtype.names:
            name=n('s2n_r')
        elif n('s2n') in data.dtype.names:
            name=n('s2n')
        elif n('s2n_w') in data.dtype.names:
            print("Using s2n_w for selections")
            name=n('s2n_w')
        else:
            name=None

        return name

    def _get_T(self, n, data, w):

        d = None
        for ntry in ['T_r','T','pars']:
            name = n(ntry)
            if name in data.dtype.names:

                if ntry=='pars':
                    d = data[name][w,4]
                else:
                    d = data[name][w]
                break

        return d

    def _get_psf_T(self, data, w):
        ns=['mcal_psfrec_T','mcal_Tpsf']

        for n in ns:
            if n in data.dtype.names:
                return data[n][w]

        return None


    def _get_s2n(self, n, data, w):
        name=self._get_s2n_name(n, data)
        if name is None:
            return None

        return data[name][w]

    def _get_size(self, n, data, w):
        name=n('pars')

        size=None
        if name in data.dtype.names:
            if data[name].shape[1] >= 6:
                size = data[name][w,4]

        return size

    def _get_flux(self, n, data, w):
        name=n('flux')

        if name in data.dtype.names:
            flux = data[name][w]
        else:
            flux=None

        return flux

    def _get_am_flux(self, n, data, w):
        name=n('am_flux')

        if name in data.dtype.names:
            flux = data[name][w]
        else:
            flux=None

        return flux

    def _get_am_flux_s2n(self, n, data, w):
        name=n('am_flux_s2n')

        if name in data.dtype.names:
            flux = data[name][w]
        else:
            flux=None

        return flux



    def _get_flux_true(self, data, w):
        if 'flux_true' in data.dtype.names:
            return data['flux_true'][w]
        else:
            return None


    def _get_flux_s2n(self, n, data, w):
        name=n('flux_s2n')
        return data[name][w]

    def _read_means(self):
        fname=self._get_means_file()
        print("reading:",fname)
        return fitsio.read(fname)

    def _write_means(self):
        fname=self._get_means_file()
        eu.ostools.makedirs_fromfile(fname)
        print("writing:",fname)
        fitsio.write(fname, self.means, extname="corr", clobber=True)
        fitsio.write(fname, self.means_nocorr, extname="nocorr")

    def _write_sums(self, run, sums):
        fname=self._get_sums_file(run)
        eu.ostools.makedirs_fromfile(fname)
        print("writing:",fname)
        fitsio.write(fname, sums, clobber=True)


    def _get_fname_extra(self, last=None, run=None, use_select_string=False):

        # name will be singular if we use args.runs
        if run is None:
            runs=self.args.runs
        else:
            runs=run

        #if len(runs) > 1:
        #    extra=runs[1:]
        #else:
        extra=[]

        if self.args.weighted:
            extra += ['weighted',self.args.weight_type]

        if self.args.preselect:
            extra += ['preselect']

        if self.args.select is not None:
            if use_select_string:
                s=self.select.replace(' ','-').replace('(','').replace(')','').replace('[','').replace(']','').replace('"','').replace("'",'')
                extra += ['select',s]
            else:
                extra += [self.args.select]

        if self.args.ntest is not None:
            extra += ['test%d' % self.args.ntest]

        if last is not None:
            extra += [last]

        if self.args.index is not None:
            extra += ['%06d' % self.args.index]

        if len(extra) > 0:
            extra = '-'.join(extra)
            extra = extra.replace('run-','')
        else:
            extra=None

        return extra



    def _get_means_file(self):

        extra=self._get_fname_extra()
        fname=files.get_means_url(self.args.runs, extra=extra)
        return fname

    def _get_sums_file(self, run, use_select_string=False):

        extra=self._get_fname_extra(
            run=run,
            use_select_string=use_select_string,
        )
        fname=files.get_sums_url(run, extra=extra)
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
            ('gsq','f8',2),
            ('wsq','f8',2),

            ('gpsf','f8',2),
            ('gpsf_orig','f8',2),

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

            ('shear_sum','f8',2),
            ('R11_sum','f8'),
            ('R22_sum','f8'),
        ]
        return dt

    def _set_select(self):
        self.select=None
        self.do_selection=False

        if self.args.select is not None:
            self.do_selection=True

            d = files.read_config(self.args.select)
            self.select = d['select'].strip()

        elif self.args.weighted:
            self.do_selection=True


    def plot_fits(self):

        means=self.means
        if means.size == 1:
            return
        else:
            import biggles
            biggles.configure('default','fontsize_min',1.5)

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

        means=self.means

        if means.size == 1:
            return
        else:
            import biggles

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

    wts=numpy.zeros(g_cov.shape[0])

    g00=g_cov[:,0,0]
    g11=g_cov[:,1,1]

    w,=numpy.where( numpy.isfinite(g00) & numpy.isfinite(g11) )
    wts[w] = 1.0/(2*args.shapenoise**2 + g00[w] + g11[w] )

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

def get_m_c_oneshear(data, nsig=2.0):

    shmeas=data['shear'][0]
    shmeas_err  = data['shear_err'][0]
    shtrue = data['shear_true'][0]

    if shtrue[1] == 0.0:
        mel=0
        cel=1
    else:
        mel=1
        cel=0

    print("errors are %g sigma" % nsig)

    fits=numpy.zeros(1, dtype=[('m','f8'),
                               ('merr','f8'),
                               ('c','f8'),
                               ('cerr','f8')])
    m=shmeas[mel]/shtrue[mel]-1.0
    merr=shmeas_err[mel]/shtrue[mel]

    c=shmeas[cel]
    cerr=shmeas_err[cel]

    fits['m']=m
    fits['merr']=merr

    fits['c']=c
    fits['cerr']=cerr


    print("m: %.3e +/- %.3e  c: %.3e +/- %.3e" % (m,nsig*merr,c,nsig*cerr))

    return fits


def add_sums(sums_in, new_sums):

    sums=sums_in.copy()
    for n in sums.dtype.names:
        if n in ['shear_true']:
            continue
        else:
            sums[n] += new_sums[n]

    return sums

def fit_m_c(data, doprint=True, onem=False, max_shear=None, nocorr_select=False, nsig=2.0):

    strue = data['shear_true']

    sdiff = data['shear'] - data['shear_true']

    serr  = data['shear_err']

    if max_shear is not None:
        stot_true = sqrt(strue[:,0]**2 + strue[:,1]**2)
        w,=where(stot_true < max_shear)
        if w.size == 0:
            raise ValueError("no shears less than %g" % max_shear)
        print("kept %d/%d with shear < %g" % (w.size,data.size,max_shear))
        strue=strue[w,:]
        sdiff=sdiff[w,:]
        serr=serr[w,:]


    m = numpy.zeros(2)
    merr = numpy.zeros(2)
    c = numpy.zeros(2)
    cerr = numpy.zeros(2)

    print("errors are %g sigma" % nsig)
    if onem:
        fits=numpy.zeros(1, dtype=[('m','f8'),
                                   ('merr','f8'),
                                   ('c1','f8'),
                                   ('c1err','f8'),
                                   ('c2','f8'),
                                   ('c2err','f8')])


        fitter=MCFitter(strue, sdiff, serr)
        fitter.dofit()
        res=fitter.get_result()

        pars=res['pars']
        perr=res['perr']
        fits['m'] = pars[0]
        fits['c1'] = pars[1]
        fits['c2'] = pars[2]
        fits['merr'] = perr[0]
        fits['c1err'] = perr[1]
        fits['c2err'] = perr[2]

        if doprint:
            print('  m:  %.3e +/- %.3e' % (pars[0],nsig*perr[0]))
            print('  c1: %.3e +/- %.3e' % (pars[1],nsig*perr[1]))
            print('  c2: %.3e +/- %.3e' % (pars[2],nsig*perr[2]))
        return fits

    fits=numpy.zeros(1, dtype=[('m','f8',2),
                               ('merr','f8',2),
                               ('c','f8',2),
                               ('cerr','f8',2),
                               ('r','f8',2)])

    for i in [0,1]:

        #c, c_err, m, m_err, covar = fitline(strue[:,i], sdiff[:,i])
        res = fitline(strue[:,i], sdiff[:,i])
        r = res['cov']/numpy.sqrt(res['slope_err']**2 * res['offset_err']**2)
        fits['m'][0,i] = res['slope']
        fits['merr'][0,i] = res['slope_err']
        fits['c'][0,i] = res['offset']
        fits['cerr'][0,i] = res['offset_err']
        fits['r'][0,i] = r

        if doprint:
            print_m_c(i+1,
                      res['slope'],
                      nsig*res['slope_err'],
                      res['offset'],
                      nsig*res['offset_err'],
                      r=r)

    return fits


