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

try:
    import reredux
except ImportError:
    pass

class SummerDT(dict):
    def __init__(self, conf, args):
        self.update(conf)
        self.args=args

        self.chunksize=1000000

        sconf=self['simc']

        self['nshear']=len(self['simc']['shear']['shears'])

        shears_array=array(self['simc']['shear']['shears'])

        # fake shear for selection effects
        sh = sqrt(shears_array[:,0]**2 + shears_array[:,1]**2)
        # in each component
        self['fake_shear']=sh.mean()/sqrt(2)

        #self['fake_shear']=0.045
        #self['fake_shear']=0.01

        print("fake shear:",self['fake_shear'])

        self._set_select()

    def go(self):

        args=self.args
        if args.fit_only:
            self.means=self._read_means()
        else:

            extra=self._get_fname_extra()

            sums, sums_select = self.do_sums()


            for i,tsums in enumerate([sums,sums_select]):
                if tsums is None:
                    continue

                if i == 0:
                    n=BNamer()
                else:
                    n=BNamer('sel')

                # averated over all shear fields
                wtot=tsums['wsum'].sum()
                self[n('R')] = tsums['R'].sum(axis=0)/wtot
                self[n('Rpsf')] = tsums['Rpsf'].sum(axis=0)/wtot

                if not args.nocorr:
                    Rnoise, Rnoise_err, Rnoise_psf, Rnoise_psf_err = \
                            get_Rnoise_line(self, args, None, extra, sums=tsums)

                    print("R:",self[n('R')])
                    print("Rnoise:",Rnoise)
                    self[n('R')] -= Rnoise
                    self[n('Rpsf')] -= Rnoise_psf
                    self[n('Rnoise')] = Rnoise
                    self[n('Rnoise_psf')] = Rnoise_psf

                self[n('Rinv')] = numpy.linalg.inv( self[n('R')] )

            self.sel=ones(2)
            n=BNamer()
            if sums_select is not None:

                n=BNamer('sel')
                sums_no_select=sums
                sums=sums_select

                if not args.nocorr_select:
                    self.sel = self.get_selection_effect()

            print("sel:",self.sel)

            g = sums['g'].copy()
            gpsf = sums['gpsf'].copy()

            winv = 1.0/sums['wsum']
            g[:,0]    *= winv
            g[:,1]    *= winv
            gpsf[:,0] *= winv
            gpsf[:,1] *= winv

            # using mean responses
            Rpsf = self[n('Rpsf')]
            Rinv = self[n('Rinv')]

            shears=self['simc']['shear']['shears']
            means=get_mean_struct(self['nshear'])
            for i in xrange(self['nshear']):

                shear_true = shears[i]

                psf_corr  = Rpsf*gpsf[i]

                gmean     = g[i]
                shear     = numpy.dot(Rinv, gmean-psf_corr)

                means['shear'][i] = shear
                means['shear_err'][i] = 1.0
                means['shear_true'][i] = shear_true

            means['shear'][:,0] *= self.sel[0]
            means['shear'][:,1] *= self.sel[1]

            self.means=means
            self._write_means()

        if args.boot:
            self.fits, self.fitsone = reredux.averaging.fit_m_c_boot(self.means)
        else:
            self.fits=reredux.averaging.fit_m_c(self.means, max_shear=args.max_shear)
            self.fitsone=reredux.averaging.fit_m_c(self.means,onem=True, max_shear=args.max_shear)

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

        if self.select is not None:
            s=self.select.replace(' ','-').replace('(','').replace(')','').replace('[','').replace(']','').replace('"','').replace("'",'')
            extra += ['select-']

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



    def get_selection_effect(self):

        sums, sums_select = self.do_sums(fake_shear=self['fake_shear'])

        sheard={}

        for i,tsums in enumerate([sums, sums_select]):
            if i==0:
                n=BNamer()
            else:
                n=BNamer('sel')

            # we added shear to all fields, so we sum over shear_index
            wtot = tsums['wsum'].sum()
            g    = tsums['g'].sum(axis=0)/wtot
            gpsf = tsums['gpsf'].sum(axis=0)/wtot

            #psf_corr  = self[n('Rpsf')]*gpsf
            #sheard[n('shear')] = numpy.dot(self[n('Rinv')], g-psf_corr)
            sheard[n('shear')] = numpy.dot(self[n('Rinv')], g)

        print("shear:",sheard['shear'])
        print("shear_sel:",sheard['shear_sel'])
        sel = sheard['shear']/sheard['shear_sel']

        return sel

    def get_selection_effect_sens(self):
        """
        calculate a full sensitivity ratio
        """
        sumsp, sums_selectp = self.do_sums(fake_shear= self['fake_shear'])
        sumsm, sums_selectm = self.do_sums(fake_shear=-self['fake_shear'])

        sheard={}

        slist=[sumsp,sumsm,sums_selectp, sums_selectm]
        nlist=['p','m','selp','selm']

        for i in xrange(len(slist)):
            tsums = slist[i]
            name=nlist[i]

            if 'sel' in name:
                n=BNamer('sel')
            else:
                n=BNamer()


            # we added shear to all fields, so we sum over shear_index
            wtot = tsums['wsum'].sum()
            g    = tsums['g'].sum(axis=0)/wtot
            gpsf = tsums['gpsf'].sum(axis=0)/wtot

            psf_corr  = self[n('Rpsf')]*gpsf

            sheard[name] = numpy.dot(self[n('Rinv')], g-psf_corr)

        print(sheard)
        Rall = sheard['p'] - sheard['m']
        Rsel = sheard['selp'] - sheard['selm']
        print("Rall:",Rall)
        print("Rsel:",Rsel)

        sel=Rall/Rsel

        return sel


    def do_sums(self, fake_shear=None):

        nrand=1
        if (fake_shear is not None
                and self.args.nrand is not None):
            nrand=self.args.nrand

        chunksize=self.chunksize
        args=self.args

        sconf=self['simc']

        if 'Rnoise_sel' in self:
            Rnoise=self['Rnoise']
            Rnoise_sel=self['Rnoise_sel']
        else:
            Rnoise=None
            Rnoise_sel=None

        sums=None
        sums_select=None
        ntot=0
        for run in args.runs:
            if args.ntest is not None and ntot > args.ntest:
                break

            fname = nsim.files.get_output_url(run, 0, 0)
            if args.d is not None:
                fname=os.path.join(args.d,os.path.basename(fname))

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
                    ntot += data.size

                    g1all=data['mcal_g'][:,0]
                    g2all=data['mcal_g'][:,1]

                    for irand in xrange(nrand):
                        if nrand > 1: print("        irand: %d/%d" % (irand+1,nrand))

                        if fake_shear is not None:
                            if self.args.nrand is None:
                                g1send,g2send=self.add_fake_shear(fake_shear,
                                                                  data,
                                                                  Rnoise=Rnoise)
                            else:
                                g1send, g2send=self.add_fake_shear_rand(fake_shear, data, Rnoise=Rnoise)
                        else:
                            g1send=g1all
                            g2send=g2all

                        sums=self.do_sums1(data, g1send, g2send, sums=sums)

                        if args.weights is not None:
                            if irand == 0:
                                if args.weights=='s2n':
                                    wts = get_s2n_weights(data['s2n_r'], args)
                                elif args.weights=='noise':
                                    wts = get_noise_weights(data['mcal_g_cov'], args)
                                else:
                                    raise ValueError("bad weight type: '%s'" % args.weights)

                            if fake_shear is not None:
                                if self.args.nrand is None:
                                    g1send,g2send=self.add_fake_shear(fake_shear,
                                                                      data,
                                                                      Rnoise=Rnoise_sel)
                                else:
                                    g1send,g2send=self.add_fake_shear_rand(fake_shear,
                                                                           data,
                                                                           Rnoise=Rnoise_sel)
                            else:
                                g1send=g1all
                                g2send=g2all

                            sums_select=self.do_sums1(data, g1send, g2send,
                                                      sums=sums_select,
                                                      weights=wts)

                        elif self.select is not None:

                            if irand == 0:
                                # only need to select once; the randomize is of the
                                # orientations not the objects used
                                logic=eval(self.select)
                                w,=numpy.where(logic)
                                print("        keeping %d/%d from cuts" % (w.size,data.size))
                                sdata=data[w]

                                g1send=sdata['mcal_g'][:,0]
                                g2send=sdata['mcal_g'][:,1]

                            # Rnoise is different, so we add shear separately for
                            # the selected data

                            if fake_shear is not None:
                                if self.args.nrand is None:
                                    g1send,g2send=self.add_fake_shear(fake_shear,
                                                                      sdata,
                                                                      Rnoise=Rnoise_sel)
                                else:
                                    g1send,g2send=self.add_fake_shear_rand(fake_shear,
                                                                           sdata,
                                                                           Rnoise=Rnoise_sel)

                            sums_select=self.do_sums1(sdata, g1send, g2send, sums=sums_select)

                    beg = beg + chunksize

                    if args.ntest is not None and ntot > args.ntest:
                        break

        if fake_shear is not None:
            # the shear is being simulated
            # we don't add the Rpsf when adding fake shear; should we?
            for n in ['Rpsf','Rdt_psf']:
                if n in sums.dtype.names:
                    sums[n] = 0.0
                    sums_select[n] = 0.0

        return sums, sums_select

    def do_sums1(self, data, g1, g2, sums=None, weights=None):
        """
        just a binner and summer, no logic here
        """
        if weights is not None:
            return self.do_sums1_weights(data, g1, g2, weights, sums=sums)

        nshear=self['nshear']
        args=self.args

        n_detrend = data['mcal_dt_Rnoise'].shape[1]

        h,rev = eu.stat.histogram(data['shear_index'],
                                  min=0,
                                  max=nshear-1,
                                  rev=True)
        nind = h.size
        assert nshear==nind

        if sums is None:
            sums=self._get_sums_struct(n_detrend)

        for i in xrange(nshear):
            if rev[i] != rev[i+1]:
                w=rev[ rev[i]:rev[i+1] ]

                t=data[w]


                sums['wsum'][i] += w.size
                sums['g'][i,0] += g1[w].sum()
                sums['g'][i,1] += g2[w].sum()
                sums['gpsf'][i] += t['mcal_gpsf'].sum(axis=0)

                sums['R'][i] += t['mcal_R'].sum(axis=0)
                sums['Rpsf'][i] += t['mcal_Rpsf'].sum(axis=0)

                if not args.nocorr:

                    sums['Rdt'][i] += t['mcal_dt_Rnoise'].sum(axis=0)

                    if not args.no_Rpsf_noise:
                        sums['Rdt_psf'][i] += t['mcal_dt_Rnoise_psf'].sum(axis=0)

        return sums

    def do_sums1_weights(self, data, g1, g2, weights, sums=None):
        """
        just a binner and summer, no logic here
        """
        nshear=self['nshear']
        args=self.args

        n_detrend = data['mcal_dt_Rnoise'].shape[1]

        h,rev = eu.stat.histogram(data['shear_index'],
                                  min=0,
                                  max=nshear-1,
                                  rev=True)
        nind = h.size
        assert nshear==nind

        if sums is None:
            sums=self._get_sums_struct(n_detrend)

        for i in xrange(nshear):
            if rev[i] != rev[i+1]:
                w=rev[ rev[i]:rev[i+1] ]

                t=data[w]

                wts=weights[w]
                wa = wts[:,newaxis]
                wa2 = wts[:,newaxis,newaxis]
                wa3 = wts[:,newaxis,newaxis,newaxis]
                wsum = wts.sum()

                sums['wsum'][i] += wsum
                #sums['g'][i] += (t['mcal_g']*wa).sum(axis=0)
                sums['g'][i,0] += (g1[w]*wts).sum()
                sums['g'][i,1] += (g2[w]*wts).sum()
                sums['gpsf'][i] += (t['mcal_gpsf']*wa).sum(axis=0)

                sums['R'][i] += (t['mcal_R']*wa2).sum(axis=0)
                sums['Rpsf'][i] += (t['mcal_Rpsf']*wa).sum(axis=0)

                if not args.nocorr:

                    sums['Rdt'][i] += (t['mcal_dt_Rnoise']*wa3).sum(axis=0)

                    if not args.no_Rpsf_noise:
                        sums['Rdt_psf'][i] += (t['mcal_dt_Rnoise_psf']*wa2).sum(axis=0)

        return sums

    def _set_select(self):
        self.select=None
        if self.args.select is not None:
            self.select = self.args.select
        elif self.args.select_from is not None:
            with open(self.args.select_from) as fobj:
                d=yaml.load(fobj)

            self.select = d['select'].strip()


    def _get_sums_struct(self, n_detrend):
        dt=self._get_sums_dt(n_detrend)
        return numpy.zeros(self['nshear'], dtype=dt)

    def _get_sums_dt(self, n_detrend):
        dt=[
            ('wsum','f8'),
            ('g','f8',2),
            ('gpsf','f8',2),
            ('R','f8',(2,2)),
            ('Rpsf','f8',2),
            ('Rdt','f8',(n_detrend,2,2)),
            ('Rdt_psf','f8',(n_detrend,2)),
        ]
        return dt

    def add_fake_shear(self, shear, data_in, Rnoise=None):
        """
        add shear using the individual R values
        """

        data=data_in.copy()

        g = data['mcal_g']

        n=g.shape[0]

        R = data['mcal_R']

        if Rnoise is not None:
            R11 = R[:,0,0].copy()
            R22 = R[:,1,1].copy()
            R11 -= Rnoise[0,0]
            R22 -= Rnoise[1,1]
        else:
            R11 = R[:,0,0]
            R22 = R[:,1,1]

        if self.args.reflect:
            tR11 = R11
            tR22 = R22
            R11 = zeros(n*2)
            R22 = zeros(n*2)
            R11[0:n] = tR11
            R11[n:]  = tR11
            R22[0:n] = tR22
            R22[n:]  = tR22

        # currently subtracting mean to remove psf anisotropy, and
        # we set Rpsf=0
        if self.args.reflect:
            g1sub = g[:,0] - g[:,0].mean()
            g2sub = g[:,1] - g[:,1].mean()

            g1 = zeros(n*2)
            g2 = zeros(n*2)
            g1[0:n] = g1sub
            g1[n:] = -g1sub
            g2[0:n] = g2sub
            g2[n:] = -g2sub

        else:
            g1 = g[:,0] - g[:,0].mean()
            g2 = g[:,1] - g[:,1].mean()

        shear1 = shear*R11
        sg1,junk = ngmix.shape.shear_reduced(g1,
                                             g2,
                                             shear1,
                                             0.0)
        shear2=shear*R22
        junk,sg2 = ngmix.shape.shear_reduced(g1,
                                             g2,
                                             0.0,
                                             shear2)

        if self.args.reflect:
            tdata=data
            data=zeros(tdata.size*2, dtype=tdata.dtype)
            data[0:n] = tdata
            data[n:] = tdata

        return sg1, sg2

    def add_fake_shear_rand(self, shear, data, Rnoise=None):
        """
        randomize the orientations

        add shear using the individual R values
        """

        g = data['mcal_g'].copy()
        n=g.shape[0]

        g1 = g[:,0] - g[:,0].mean()
        g2 = g[:,1] - g[:,1].mean()

        gtot = sqrt(g1**2 + g2**2)
        twotheta = numpy.random.uniform(low=0.0, high=2.0*numpy.pi, size=n)

        numpy.cos(twotheta, g1)
        numpy.sin(twotheta, g2)

        g1 *= gtot
        g2 *= gtot

        R = data['mcal_R']

        if Rnoise is not None:
            R11 = R[:,0,0].copy()
            R22 = R[:,1,1].copy()
            R11 -= Rnoise[0,0]
            R22 -= Rnoise[1,1]
        else:
            R11 = R[:,0,0]
            R22 = R[:,1,1]

        shear1 = shear*R11
        sg1,junk = ngmix.shape.shear_reduced(g1,
                                             g2,
                                             shear1,
                                             0.0)
        shear2=shear*R22
        junk,sg2 = ngmix.shape.shear_reduced(g1,
                                             g2,
                                             0.0,
                                             shear2)

        return sg1, sg2

    def add_fake_shear_noR(self, data_in):
        """
        add shear without response
        """

        data=data_in.copy()

        shear=self['fake_shear']

        g = data['mcal_g']

        g1 = g[:,0] - g[:,0].mean()
        g2 = g[:,1] - g[:,1].mean()

        print("        making sheared")
        sg1,junk = ngmix.shape.shear_reduced(g1,
                                             g2,
                                             shear,
                                             0.0)
        junk,sg2 = ngmix.shape.shear_reduced(g1,
                                             g2,
                                             0.0,
                                             shear)

        data['mcal_g'][:,0] = sg1
        data['mcal_g'][:,1] = sg2

        return data

    def plot_fits(self):
        import biggles

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
                xlabel='shear%d true' % (i+1,),
                ylabel='shear%d diff' % (i+1,),
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

class Summer(SummerDT):
    """
    metacal summer without detrend, either plain metacal
    or with fixnoise
    """
    def go(self):

        if self.args.fit_only:
            self.means=self._read_means()
        else:

            extra=self._get_fname_extra()

            sums, sums_select = self.do_sums()

            args=self.args

            for i,tsums in enumerate([sums,sums_select]):
                if tsums is None:
                    continue

                if i == 0:
                    n=BNamer()
                else:
                    n=BNamer('sel')

                # averated over all shear fields
                wtot=tsums['wsum'].sum()
                self[n('R')] = tsums['R'].sum(axis=0)/wtot
                self[n('Rpsf')] = tsums['Rpsf'].sum(axis=0)/wtot
                self[n('Rinv')] = numpy.linalg.inv( self[n('R')] )
                print("R:",self[n('R')])

            self.sel=ones(2)
            n=BNamer()
            if sums_select is not None:

                n=BNamer('sel')
                sums_no_select=sums
                sums=sums_select

                if not args.nocorr_select:
                    self.sel = self.get_selection_effect()

            print("sel:",self.sel)

            g = sums['g'].copy()
            gpsf = sums['gpsf'].copy()

            winv = 1.0/sums['wsum']
            g[:,0]    *= winv
            g[:,1]    *= winv
            gpsf[:,0] *= winv
            gpsf[:,1] *= winv

            # using mean responses
            Rpsf = self[n('Rpsf')]
            Rinv = self[n('Rinv')]

            shears=self['simc']['shear']['shears']
            means=get_mean_struct(self['nshear'])
            for i in xrange(self['nshear']):

                shear_true = shears[i]

                psf_corr  = Rpsf*gpsf[i]

                gmean = g[i]
                shear = numpy.dot(Rinv, gmean-psf_corr)

                means['shear'][i] = shear
                means['shear_err'][i] = 1.0
                means['shear_true'][i] = shear_true

            means['shear'][:,0] *= self.sel[0]
            means['shear'][:,1] *= self.sel[1]

            self.means=means
            self._write_means()

        if self.args.boot:
            self.fits, self.fitsone = reredux.averaging.fit_m_c_boot(self.means)
        else:
            self.fits=reredux.averaging.fit_m_c(self.means)
            self.fitsone=reredux.averaging.fit_m_c(self.means,onem=True)

    def do_sums(self, fake_shear=None):

        chunksize=self.chunksize
        args=self.args

        sconf=self['simc']

        sums=None
        sums_select=None
        ntot=0
        for run in args.runs:
            if args.ntest is not None and ntot > args.ntest:
                break

            fname = nsim.files.get_output_url(run, 0, 0)
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
                    ntot += data.size

                    g1all=data['mcal_g'][:,0]
                    g2all=data['mcal_g'][:,1]

                    if fake_shear is not None:
                        print("    running add fake shear")
                        g1send, g2send=self.add_fake_shear(fake_shear, data)
                    else:
                        g1send=g1all
                        g2send=g2all

                    sums=self.do_sums1(data, g1send, g2send, sums=sums)

                    if args.weights is not None:
                        if args.weights=='s2n':
                            wts = get_s2n_weights(data['s2n_r'], args)
                        elif args.weights=='noise':
                            wts = get_noise_weights(data['g_cov'], args)
                        else:
                            raise ValueError("bad weight type: '%s'" % args.weights)

                        if fake_shear is not None:
                            g1send,g2send=self.add_fake_shear(fake_shear, data)
                        else:
                            g1send=g1all
                            g2send=g2all

                        sums_select=self.do_sums1(data, g1send, g2send,
                                                  sums=sums_select,
                                                  weights=wts)

                    elif self.select is not None:

                        logic=eval(self.select)
                        w,=numpy.where(logic)
                        print("        keeping %d/%d from cuts" % (w.size,data.size))
                        sdata=data[w]

                        g1send=sdata['mcal_g'][:,0]
                        g2send=sdata['mcal_g'][:,1]

                        # Rnoise is different, so we add shear separately for
                        # the selected data

                        if fake_shear is not None:
                            print("    running add fake shear selected")
                            g1send,g2send=self.add_fake_shear(fake_shear, sdata)

                        sums_select=self.do_sums1(sdata, g1send, g2send, sums=sums_select)

                    beg = beg + chunksize

                    if args.ntest is not None and ntot > args.ntest:
                        break

        if fake_shear is not None:
            # the shear is being simulated
            # we don't add the Rpsf when adding fake shear
            for n in ['Rpsf']:
                sums[n] = 0.0
                sums_select[n] = 0.0


        return sums, sums_select

    def do_sums1(self, data, g1, g2, sums=None, weights=None):
        """
        just a binner and summer, no logic here
        """
        if weights is not None:
            return self.do_sums1_weights(data, g1, g2, weights, sums=sums)

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

        for i in xrange(nshear):
            if rev[i] != rev[i+1]:
                w=rev[ rev[i]:rev[i+1] ]

                t=data[w]


                sums['wsum'][i] += w.size
                sums['g'][i,0] += g1[w].sum()
                sums['g'][i,1] += g2[w].sum()
                sums['gpsf'][i] += t['mcal_gpsf'].sum(axis=0)

                sums['R'][i] += t['mcal_R'].sum(axis=0)
                sums['Rpsf'][i] += t['mcal_Rpsf'].sum(axis=0)

        return sums

    def do_sums1_weights(self, data, g1, g2, weights, sums=None):
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

        for i in xrange(nshear):
            if rev[i] != rev[i+1]:
                w=rev[ rev[i]:rev[i+1] ]

                t=data[w]

                wts=weights[w]
                wa = wts[:,newaxis]
                wa2 = wts[:,newaxis,newaxis]
                wsum = wts.sum()

                sums['wsum'][i] += wsum
                #sums['g'][i] += (t['mcal_g']*wa).sum(axis=0)
                sums['g'][i,0] += (g1[w]*wts).sum()
                sums['g'][i,1] += (g2[w]*wts).sum()
                sums['gpsf'][i] += (t['mcal_gpsf']*wa).sum(axis=0)

                sums['R'][i] += (t['mcal_R']*wa2).sum(axis=0)
                sums['Rpsf'][i] += (t['mcal_Rpsf']*wa).sum(axis=0)

        return sums

    def add_fake_shear(self, shear, data):
        """
        add shear using the individual R values
        """

        print("    doing add fake shear")

        g = data['mcal_g'].copy()
        n=g.shape[0]

        g1 = g[:,0] - g[:,0].mean()
        g2 = g[:,1] - g[:,1].mean()

        R = data['mcal_R']
        R11 = R[:,0,0]
        R22 = R[:,1,1]

        shear1 = shear*R11
        sg1,junk = ngmix.shape.shear_reduced(g1,
                                             g2,
                                             shear1,
                                             0.0)
        shear2=shear*R22
        junk,sg2 = ngmix.shape.shear_reduced(g1,
                                             g2,
                                             0.0,
                                             shear2)

        return sg1, sg2



    def _get_sums_struct(self):
        dt=self._get_sums_dt()
        return numpy.zeros(self['nshear'], dtype=dt)

    def _get_sums_dt(self):
        dt=[
            ('wsum','f8'),
            ('g','f8',2),
            ('gpsf','f8',2),
            ('R','f8',(2,2)),
            ('Rpsf','f8',2),
        ]
        return dt


class SummerNocorr(SummerDT):
    """
    no R corrections at all
    """
    def go(self):

        if self.args.fit_only:
            self.means=self._read_means()
        else:

            extra=self._get_fname_extra()

            sums, sums_select = self.do_sums()

            args=self.args

            for i,tsums in enumerate([sums,sums_select]):
                if tsums is None:
                    continue

                if i == 0:
                    n=BNamer()
                else:
                    n=BNamer('sel')

                # averated over all shear fields
                wtot=tsums['wsum'].sum()

            self.sel=ones(2)
            n=BNamer()
            if sums_select is not None:

                n=BNamer('sel')
                sums_no_select=sums
                sums=sums_select

                if not args.nocorr_select:
                    self.sel = self.get_selection_effect()

            print("sel:",self.sel)

            g = sums['g'].copy()

            winv = 1.0/sums['wsum']
            g[:,0]    *= winv
            g[:,1]    *= winv

            shears=self['simc']['shear']['shears']
            means=get_mean_struct(self['nshear'])
            for i in xrange(self['nshear']):

                shear_true = shears[i]

                shear = g[i]

                means['shear'][i] = shear
                means['shear_err'][i] = 1.0
                means['shear_true'][i] = shear_true

            means['shear'][:,0] *= self.sel[0]
            means['shear'][:,1] *= self.sel[1]

            self.means=means
            self._write_means()

        if self.args.boot:
            self.fits, self.fitsone = reredux.averaging.fit_m_c_boot(self.means)
        else:
            self.fits=reredux.averaging.fit_m_c(self.means)
            self.fitsone=reredux.averaging.fit_m_c(self.means,onem=True)

    def do_sums(self, fake_shear=None):

        nrand=1
        if fake_shear is not None:
            nrand=self.args.nrand

        chunksize=self.chunksize
        args=self.args

        sconf=self['simc']

        sums=None
        sums_select=None
        ntot=0
        for run in args.runs:
            if args.ntest is not None and ntot > args.ntest:
                break

            fname = nsim.files.get_output_url(run, 0, 0)
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
                    ntot += data.size

                    g1all=data['g'][:,0]
                    g2all=data['g'][:,1]

                    for irand in xrange(nrand):
                        if nrand > 1: print("        irand: %d/%d" % (irand+1,nrand))

                        if fake_shear is not None:
                            #tmp_data=self.add_fake_shear(fake_shear, data, self['Rnoise'])
                            # use sheared values
                            g1send, g2send=self.add_fake_shear_rand(fake_shear, data)
                        else:
                            g1send=g1all
                            g2send=g2all

                        sums=self.do_sums1(data, g1send, g2send, sums=sums)

                        if args.weights is not None:
                            if irand == 0:
                                if args.weights=='s2n':
                                    wts = get_s2n_weights(data['s2n_r'], args)
                                elif args.weights=='noise':
                                    wts = get_noise_weights(data['g_cov'], args)
                                else:
                                    raise ValueError("bad weight type: '%s'" % args.weights)

                            if fake_shear is not None:
                                #tmp_data=self.add_fake_shear(fake_shear, data, self['Rnoise_sel'])
                                # use sheared values, with different Rnoise
                                g1send,g2send=self.add_fake_shear_rand(fake_shear,
                                                                       data)
                            else:
                                g1send=g1all
                                g2send=g2all

                            sums_select=self.do_sums1(data, g1send, g2send,
                                                      sums=sums_select,
                                                      weights=wts)

                        elif self.select is not None:

                            if irand == 0:
                                logic=eval(self.select)
                                w,=numpy.where(logic)
                                print("        keeping %d/%d from cuts" % (w.size,data.size))
                                sdata=data[w]

                                g1send=sdata['g'][:,0]
                                g2send=sdata['g'][:,1]

                            # Rnoise is different, so we add shear separately for
                            # the selected data

                            if fake_shear is not None:
                                #tmp_sdata=self.add_fake_shear(fake_shear, sdata, self['Rnoise_sel'])
                                g1send,g2send=self.add_fake_shear_rand(fake_shear,
                                                                       sdata)

                            sums_select=self.do_sums1(sdata, g1send, g2send, sums=sums_select)

                    beg = beg + chunksize

                    if args.ntest is not None and ntot > args.ntest:
                        break

        return sums, sums_select

    def do_sums1(self, data, g1, g2, sums=None, weights=None):
        """
        just a binner and summer, no logic here
        """
        if weights is not None:
            return self.do_sums1_weights(data, g1, g2, weights, sums=sums)

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

        for i in xrange(nshear):
            if rev[i] != rev[i+1]:
                w=rev[ rev[i]:rev[i+1] ]

                t=data[w]


                sums['wsum'][i] += w.size
                #sums['g'][i] += t['mcal_g'].sum(axis=0)
                sums['g'][i,0] += g1[w].sum()
                sums['g'][i,1] += g2[w].sum()

        return sums

    def _get_sums_struct(self):
        dt=self._get_sums_dt()
        return numpy.zeros(self['nshear'], dtype=dt)

    def _get_sums_dt(self):
        dt=[
            ('wsum','f8'),
            ('g','f8',2),
        ]
        return dt


class SummerPCal(dict):
    def __init__(self, conf, args):
        self.update(conf)
        self.args=args

        self.chunksize=1000000

        sconf=self['simc']

        self['nshear']=len(self['simc']['shear']['shears'])

        shears_array=array(self['simc']['shear']['shears'])

        '''
        # fake shear for selection effects
        sh = sqrt(shears_array[:,0]**2 + shears_array[:,1]**2)
        # in each component
        self['fake_shear']=sh.mean()/sqrt(2)

        #self['fake_shear']=0.045
        '''
        self['fake_shear']=0.01

        print("fake shear:",self['fake_shear'])

        self._set_select()

    def go(self):

        if self.args.fit_only:
            self.means=self._read_means()
        else:

            extra=self._get_fname_extra()

            sums, sums_select = self.do_sums()

            args=self.args

            for i,tsums in enumerate([sums,sums_select]):
                if tsums is None:
                    continue

                if i == 0:
                    n=BNamer()
                else:
                    n=BNamer('sel')

                # averated over all shear fields
                wtot=tsums['wsum'].sum()
                self[n('R')] = tsums['R'].sum(axis=0)/wtot
                self[n('Rpsf')] = tsums['Rpsf'].sum(axis=0)/wtot

                if not args.nocorr:
                    Rnoise, Rnoise_err, Rnoise_psf, Rnoise_psf_err = \
                            get_Rnoise_line(self, args, None, extra, sums=tsums)

                    print("R:",self[n('R')])
                    print("Rnoise:",Rnoise)
                    self[n('R')] -= Rnoise
                    self[n('Rpsf')] -= Rnoise_psf
                    self[n('Rnoise')] = Rnoise
                    self[n('Rnoise_psf')] = Rnoise_psf

                self[n('Rinv')] = numpy.linalg.inv( self[n('R')] )

            self.sel=ones(2)
            n=BNamer()
            if sums_select is not None:

                n=BNamer('sel')
                sums_no_select=sums
                sums=sums_select

                if not args.nocorr_select:
                    self.sel = self.get_selection_effect()

            print("sel:",self.sel)

            g = sums['g'].copy()
            gpsf = sums['gpsf'].copy()

            winv = 1.0/sums['wsum']
            g[:,0]    *= winv
            g[:,1]    *= winv
            gpsf[:,0] *= winv
            gpsf[:,1] *= winv

            # using mean responses
            Rpsf = self[n('Rpsf')]
            Rinv = self[n('Rinv')]

            shears=self['simc']['shear']['shears']
            means=get_mean_struct(self['nshear'])
            for i in xrange(self['nshear']):

                shear_true = shears[i]

                psf_corr  = Rpsf*gpsf[i]

                gmean     = g[i]
                shear     = numpy.dot(Rinv, gmean-psf_corr)

                means['shear'][i] = shear
                means['shear_err'][i] = 1.0
                means['shear_true'][i] = shear_true

            means['shear'][:,0] *= self.sel[0]
            means['shear'][:,1] *= self.sel[1]

            self.means=means
            self._write_means()

        if self.args.boot:
            self.fits, self.fitsone = reredux.averaging.fit_m_c_boot(self.means)
        else:
            self.fits=reredux.averaging.fit_m_c(self.means)
            self.fitsone=reredux.averaging.fit_m_c(self.means,onem=True)

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

        if self.select is not None:
            s=self.select.replace(' ','-').replace('(','').replace(')','').replace('[','').replace(']','').replace('"','').replace("'",'')
            extra += ['select-']

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



    def get_selection_effect(self):

        sums, sums_select = self.do_sums(fake_shear=self['fake_shear'])

        sheard={}

        for i,tsums in enumerate([sums, sums_select]):
            if i==0:
                n=BNamer()
            else:
                n=BNamer('sel')

            # we added shear to all fields, so we sum over shear_index
            wtot = tsums['wsum'].sum()
            g    = tsums['g'].sum(axis=0)/wtot
            gpsf = tsums['gpsf'].sum(axis=0)/wtot

            #psf_corr  = self[n('Rpsf')]*gpsf
            #sheard[n('shear')] = numpy.dot(self[n('Rinv')], g-psf_corr)
            sheard[n('shear')] = numpy.dot(self[n('Rinv')], g)

        print("shear:",sheard['shear'])
        print("shear_sel:",sheard['shear_sel'])
        sel = sheard['shear']/sheard['shear_sel']

        return sel

    def get_selection_effect_sens(self):
        """
        calculate a full sensitivity ratio
        """
        sumsp, sums_selectp = self.do_sums(fake_shear= self['fake_shear'])
        sumsm, sums_selectm = self.do_sums(fake_shear=-self['fake_shear'])

        sheard={}

        slist=[sumsp,sumsm,sums_selectp, sums_selectm]
        nlist=['p','m','selp','selm']

        for i in xrange(len(slist)):
            tsums = slist[i]
            name=nlist[i]

            if 'sel' in name:
                n=BNamer('sel')
            else:
                n=BNamer()


            # we added shear to all fields, so we sum over shear_index
            wtot = tsums['wsum'].sum()
            g    = tsums['g'].sum(axis=0)/wtot
            gpsf = tsums['gpsf'].sum(axis=0)/wtot

            psf_corr  = self[n('Rpsf')]*gpsf

            sheard[name] = numpy.dot(self[n('Rinv')], g-psf_corr)

        print(sheard)
        Rall = sheard['p'] - sheard['m']
        Rsel = sheard['selp'] - sheard['selm']
        print("Rall:",Rall)
        print("Rsel:",Rsel)

        sel=Rall/Rsel

        return sel


    def do_sums(self, fake_shear=None):

        nrand=1
        if fake_shear is not None:
            nrand=self.args.nrand

        chunksize=self.chunksize
        args=self.args

        sconf=self['simc']

        sums=None
        sums_select=None
        ntot=0
        for run in args.runs:
            if args.ntest is not None and ntot > args.ntest:
                break

            fname = nsim.files.get_output_url(run, 0, 0)
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
                    ntot += data.size

                    g1all=data['mcal_g'][:,0]
                    g2all=data['mcal_g'][:,1]

                    for irand in xrange(nrand):
                        if nrand > 1: print("        irand: %d/%d" % (irand+1,nrand))

                        if fake_shear is not None:
                            #tmp_data=self.add_fake_shear(fake_shear, data, self['Rnoise'])
                            # use sheared values
                            g1send, g2send=self.add_fake_shear_rand(fake_shear, data, self['Rnoise'])
                        else:
                            g1send=g1all
                            g2send=g2all

                        sums=self.do_sums1(data, g1send, g2send, sums=sums)

                        if args.weights is not None:
                            if nrand == 0:
                                if args.weights=='s2n':
                                    wts = get_s2n_weights(tmp_data['s2n_r'], args)
                                elif args.weights=='noise':
                                    wts = get_noise_weights(tmp_data['mcal_g_cov'], args)
                                else:
                                    raise ValueError("bad weight type: '%s'" % args.weights)

                            if fake_shear is not None:
                                #tmp_data=self.add_fake_shear(fake_shear, data, self['Rnoise_sel'])
                                # use sheared values, with different Rnoise
                                g1send,g2send=self.add_fake_shear_rand(fake_shear,
                                                                       data,
                                                                       self['Rnoise_sel'])
                            else:
                                g1send=g1all
                                g2send=g2all

                            sums_select=self.do_sums1(data, g1send, g2send,
                                                      sums=sums_select,
                                                      weights=wts)

                        elif self.select is not None:

                            if irand == 0:
                                logic=eval(self.select)
                                w,=numpy.where(logic)
                                print("        keeping %d/%d from cuts" % (w.size,data.size))
                                sdata=data[w]

                                g1send=sdata['mcal_g'][:,0]
                                g2send=sdata['mcal_g'][:,1]

                            # Rnoise is different, so we add shear separately for
                            # the selected data

                            if fake_shear is not None:
                                #tmp_sdata=self.add_fake_shear(fake_shear, sdata, self['Rnoise_sel'])
                                g1send,g2send=self.add_fake_shear_rand(fake_shear,
                                                                       sdata,
                                                                       self['Rnoise_sel'])

                            sums_select=self.do_sums1(sdata, g1send, g2send, sums=sums_select)

                    beg = beg + chunksize

                    if args.ntest is not None and ntot > args.ntest:
                        break

        if fake_shear is not None:
            # the shear is being simulated
            # we don't add the Rpsf when adding fake shear; should we?
            for n in ['Rpsf','Rdt_psf']:
                sums[n] = 0.0
                sums_select[n] = 0.0

        return sums, sums_select

    def do_sums1(self, data, g1, g2, sums=None, weights=None):
        """
        just a binner and summer, no logic here
        """
        if weights is not None:
            return self.do_sums1_weights(data, g1, g2, weights, sums=sums)

        nshear=self['nshear']
        args=self.args

        n_detrend = data['mcal_dt_Rnoise'].shape[1]

        h,rev = eu.stat.histogram(data['shear_index'],
                                  min=0,
                                  max=nshear-1,
                                  rev=True)
        nind = h.size
        assert nshear==nind

        if sums is None:
            sums=self._get_sums_struct(n_detrend)

        for i in xrange(nshear):
            if rev[i] != rev[i+1]:
                w=rev[ rev[i]:rev[i+1] ]

                t=data[w]


                sums['wsum'][i] += w.size
                #sums['g'][i] += t['mcal_g'].sum(axis=0)
                sums['g'][i,0] += g1[w].sum()
                sums['g'][i,1] += g2[w].sum()
                sums['gpsf'][i] += t['mcal_gpsf'].sum(axis=0)

                sums['R'][i] += t['mcal_R'].sum(axis=0)
                sums['Rpsf'][i] += t['mcal_Rpsf'].sum(axis=0)

                if not args.nocorr:

                    sums['Rdt'][i] += t['mcal_dt_Rnoise'].sum(axis=0)

                    if not args.no_Rpsf_noise:
                        sums['Rdt_psf'][i] += t['mcal_dt_Rnoise_psf'].sum(axis=0)

        return sums

    def do_sums1_weights(self, data, g1, g2, weights, sums=None):
        """
        just a binner and summer, no logic here
        """
        nshear=self['nshear']
        args=self.args

        n_detrend = data['mcal_dt_Rnoise'].shape[1]

        h,rev = eu.stat.histogram(data['shear_index'],
                                  min=0,
                                  max=nshear-1,
                                  rev=True)
        nind = h.size
        assert nshear==nind

        if sums is None:
            sums=self._get_sums_struct(n_detrend)

        for i in xrange(nshear):
            if rev[i] != rev[i+1]:
                w=rev[ rev[i]:rev[i+1] ]

                t=data[w]

                wts=weights[w]
                wa = wts[:,newaxis]
                wa2 = wts[:,newaxis,newaxis]
                wa3 = wts[:,newaxis,newaxis,newaxis]
                wsum = wts.sum()

                sums['wsum'][i] += wsum
                #sums['g'][i] += (t['mcal_g']*wa).sum(axis=0)
                sums['g'][i,0] += (g1[w]*wts).sum()
                sums['g'][i,1] += (g2[w]*wts).sum()
                sums['gpsf'][i] += (t['mcal_gpsf']*wa).sum(axis=0)

                sums['R'][i] += (t['mcal_R']*wa2).sum(axis=0)
                sums['Rpsf'][i] += (t['mcal_Rpsf']*wa).sum(axis=0)

                if not args.nocorr:

                    sums['Rdt'][i] += (t['mcal_dt_Rnoise']*wa3).sum(axis=0)

                    if not args.no_Rpsf_noise:
                        sums['Rdt_psf'][i] += (t['mcal_dt_Rnoise_psf']*wa2).sum(axis=0)

        return sums

    def _set_select(self):
        self.select=None
        if self.args.select is not None:
            self.select = self.args.select
        elif self.args.select_from is not None:
            with open(self.args.select_from) as fobj:
                d=yaml.load(fobj)

            self.select = d['select'].strip()


    def _get_sums_struct(self, n_detrend):
        dt=self._get_sums_dt(n_detrend)
        return numpy.zeros(self['nshear'], dtype=dt)

    def _get_sums_dt(self, n_detrend):
        dt=[
            ('wsum','f8'),
            ('g','f8',2),
            ('gpsf','f8',2),
            ('R','f8',(2,2)),
            ('Rpsf','f8',2),
            ('Rdt','f8',(n_detrend,2,2)),
            ('Rdt_psf','f8',(n_detrend,2)),
        ]
        return dt

    def add_fake_shear(self, shear, data_in, Rnoise):
        """
        add shear using the individual R values
        """

        #print("            shear:",shear)
        data=data_in.copy()

        g = data['mcal_g']

        n=g.shape[0]

        R = data['mcal_R']
        R11 = R[:,0,0].copy()
        R22 = R[:,1,1].copy()

        R11 -= Rnoise[0,0]
        R22 -= Rnoise[1,1]

        if self.args.reflect:
            tR11 = R11
            tR22 = R22
            R11 = zeros(n*2)
            R22 = zeros(n*2)
            R11[0:n] = tR11
            R11[n:]  = tR11
            R22[0:n] = tR22
            R22[n:]  = tR22

        # currently subtracting mean to remove psf anisotropy, and
        # we set Rpsf=0
        if self.args.reflect:
            g1sub = g[:,0] - g[:,0].mean()
            g2sub = g[:,1] - g[:,1].mean()

            g1 = zeros(n*2)
            g2 = zeros(n*2)
            g1[0:n] = g1sub
            g1[n:] = -g1sub
            g2[0:n] = g2sub
            g2[n:] = -g2sub

        else:
            g1 = g[:,0] - g[:,0].mean()
            g2 = g[:,1] - g[:,1].mean()

        shear1 = shear*R11
        #print("        shear1 min max:",shear1.min(),shear1.max())
        sg1,junk = ngmix.shape.shear_reduced(g1,
                                             g2,
                                             shear1,
                                             0.0)
        shear2=shear*R22
        #print("        shear2 min max:",shear2.min(),shear2.max())
        junk,sg2 = ngmix.shape.shear_reduced(g1,
                                             g2,
                                             0.0,
                                             shear2)

        if self.args.reflect:
            tdata=data
            data=zeros(tdata.size*2, dtype=tdata.dtype)
            data[0:n] = tdata
            data[n:] = tdata

        data['mcal_g'][:,0] = sg1
        data['mcal_g'][:,1] = sg2

        return data

    def add_fake_shear_rand(self, shear, data, Rnoise):
        """
        add shear using the individual R values
        """

        g = data['mcal_g'].copy()
        n=g.shape[0]

        g1 = g[:,0] - g[:,0].mean()
        g2 = g[:,1] - g[:,1].mean()

        gtot = sqrt(g1**2 + g2**2)
        twotheta = numpy.random.uniform(low=0.0, high=2.0*numpy.pi, size=n)

        numpy.cos(twotheta, g1)
        numpy.sin(twotheta, g2)

        g1 *= gtot
        g2 *= gtot

        R = data['mcal_R']
        R11 = R[:,0,0].copy()
        R22 = R[:,1,1].copy()

        R11 -= Rnoise[0,0]
        R22 -= Rnoise[1,1]

        shear1 = shear*R11
        sg1,junk = ngmix.shape.shear_reduced(g1,
                                             g2,
                                             shear1,
                                             0.0)
        shear2=shear*R22
        junk,sg2 = ngmix.shape.shear_reduced(g1,
                                             g2,
                                             0.0,
                                             shear2)

        return sg1, sg2

    def add_fake_shear_noR(self, data_in):
        """
        add shear without response
        """

        data=data_in.copy()

        shear=self['fake_shear']

        g = data['mcal_g']

        g1 = g[:,0] - g[:,0].mean()
        g2 = g[:,1] - g[:,1].mean()

        print("        making sheared")
        sg1,junk = ngmix.shape.shear_reduced(g1,
                                             g2,
                                             shear,
                                             0.0)
        junk,sg2 = ngmix.shape.shear_reduced(g1,
                                             g2,
                                             0.0,
                                             shear)

        data['mcal_g'][:,0] = sg1
        data['mcal_g'][:,1] = sg2

        return data

    def plot_fits(self):
        import biggles

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
                xlabel='shear%d true' % (i+1,),
                ylabel='shear%d diff' % (i+1,),
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

def get_Rnoise_line(conf, args, data, extra, weights=None, sums=None):

    noise0 = conf['simc']['noise']
    factors=array(conf['detrend_factors'])
    target_noises = noise0*factors
    ndiff = target_noises - noise0
    xvals = 2*noise0*ndiff

    if not args.no_Rpsf_noise:
        do_psf=True
    else:
        do_psf=False

    if sums is not None:
        # sum over all shear fields
        winv = 1.0/sums['wsum'].sum()
        Rdt = sums['Rdt'].sum(axis=0) * winv
        Rdt_psf = sums['Rdt_psf'].sum(axis=0) * winv

    else:
        if weights is not None:
            wsum=weights.sum()

            R = data['mcal_dt_Rnoise']
            Rshape=R.shape[1:]
            Rdt = numpy.zeros(Rshape)


            for i in xrange(Rshape[0]):
                for j in xrange(Rshape[1]):
                    for k in xrange(Rshape[2]):
                        Rdt[i,j,k] = (R[:,i,j,k]*weights).sum()/wsum

            if do_psf:
                Rpsf = data['mcal_dt_Rnoise_psf']
                Rpsf_shape=Rpsf.shape[1:]

                Rdt_psf = numpy.zeros(Rpsf_shape)

                for i in xrange(Rpsf_shape[0]):
                    for j in xrange(Rpsf_shape[1]):
                        Rdt_psf[i,j] = (Rpsf[:,i,j]*weights).sum()/wsum

        else:
            Rdt = data['mcal_dt_Rnoise'].mean(axis=0)

            if do_psf:
                Rdt_psf = data['mcal_dt_Rnoise_psf'].mean(axis=0)

    A = zeros( (2,2) )
    Aerr = zeros( (2,2) )
    Apsf = zeros(2)
    Apsf_err = zeros(2)

    p='%s (%.3g +/- %.3g) + (%.3g +/ %.3g) deltan'
    for i in xrange(2):

        if do_psf:
            res = fitline(xvals, Rdt_psf[:,i])
            Apsf[i] = res['slope']
            Apsf_err[i] = res['slope_err']

            textra = 'Rnoise-detrend-Rpsf%d' % (i+1)
            if extra is not None:
                textra = '%s-%s' % (extra, textra)

            plot_line_fit(
                args,
                textra,
                xvals, Rdt_psf[:,i],res,
                r'$2 n \Delta n$',
                r'$\Delta R^{PSF}_%d$' % (i+1),
                label_error=False,
            )

        for j in xrange(2):
            res = fitline(xvals, Rdt[:,i,j])
            A[i,j] = res['slope']
            Aerr[i,j] = res['slope_err']

            n='A[%d,%d]' % (i+1,j+1)
            s=res['slope']
            serr=res['slope_err']
            o=res['offset']
            oerr=res['offset_err']
            print(p % (n,o,oerr,s,serr))

            textra='Rnoise-detrend-R%d%d' % (i+1,j+1)
            if extra is not None:
                textra = '%s-%s' % (extra, textra)
            plot_line_fit(
                args,
                textra,
                xvals, Rdt[:,i,j],res,
                r'$2 n \Delta n$',
                r'$\Delta R_{%d,%d}$' % (i+1,j+1),
                label_error=False,
            )

    Rnoise = A*noise0**2
    Rnoise_err = Aerr*noise0**2

    Rnoise_psf = Apsf*noise0**2
    Rnoise_psf_err = Apsf_err*noise0**2

    print_Rs(Rnoise, Rnoise_err, Rnoise_psf, Rnoise_psf_err, type='noise')

    return Rnoise, Rnoise_err, Rnoise_psf, Rnoise_psf_err

class BNamer(object):
    """
    create strings with a specified front prefix
    """
    def __init__(self, back=None):
        self.back=back
    def __call__(self, name):
        if self.back is None or self.back=='':
            return name
        else:
            return '%s_%s' % (name,self.back)


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


