from __future__ import print_function
import os, sys
import fitsio
import numpy
import esutil as eu
from esutil.numpy_util import between
import galsim

class DiscreteSampler(object):
    def __init__(self, vals, rng=None):
        self.vals=vals
        self.nvals=len(vals)

        if rng is None:
            rng=numpy.random.RandomState()
        self.rng=rng

    def sample(self):
        #ri=numpy.random.randint(0, self.nvals)
        ri=self.rng.randint(0, self.nvals)
        return self.vals[ri]

class DiscreteHLRFluxSampler(object):
    def __init__(self, data, rng=None):
        self.data=data
        self.ndata=len(data)

        if rng is None:
            rng=numpy.random.RandomState()
        self.rng=rng

    def sample(self):
        ri=self.rng.randint(0, self.ndata)

        r50=self.data['hlr'][ri,0]
        flux=self.data['flux'][ri,0]

        return {'r50':r50, 'flux':flux}


class PowerLaw(object):
    def __init__(self, index, xmin, xmax, npts=10000):

        func = lambda x:  x**index

        self.pdf = eu.random.Generator(
            func,
            xrange=[xmin,xmax],
            nx=npts,
        )

    def sample(self, n=None):
        return self.pdf.sample(n)

class SeparableShapeR50FluxPDF(object):
    """
    separable
    """
    def __init__(self, r50_pdf, flux_pdf, g_pdf=None):
        self.r50_pdf=r50_pdf
        self.flux_pdf=flux_pdf
        self.g_pdf=g_pdf

    def sample(self):

        if self.g_pdf is not None:
            g1,g2 = self.g_pdf.sample2d()
        else:
            g1,g2=0.0,0.0
        r50=self.r50_pdf.sample()
        flux=self.flux_pdf.sample()

        return g1, g2, r50, flux

class ShapeJointR50FluxPDF(object):
    def __init__(self, r50_flux_pdf, g_pdf=None):
        self.r50_flux_pdf=r50_flux_pdf
        self.g_pdf=g_pdf

    def sample(self):

        if self.g_pdf is not None:
            g1,g2 = self.g_pdf.sample2d()
        else:
            g1,g2=0.0,0.0
        r50,flux=self.r50_flux_pdf.sample()
        return g1, g2, r50, flux


class CosmosR50Flux(object):
    def __init__(self, r50_range, flux_range):
        self.r50_range=r50_range
        self.flux_range=flux_range

        self.r50_sanity_range=0.05,2.0
        self.flux_sanity_range=0.5,100.0
        self.kde_factor=0.01

        self._load_data()
        self._make_kde()

    def sample(self, size=None):
        """
        get [r50, flux] or [:, r50_flux]
        """
        if size is None:
            size=1
            is_scalar=True
        else:
            is_scalar=False

        r50min,r50max=self.r50_range
        fmin,fmax=self.flux_range

        data=numpy.zeros( (size,2) )

        ngood=0
        nleft=data.shape[0]
        while nleft > 0:
            r=self.kde.resample(size=nleft).T

            w,=numpy.where(
                between(r[:,0], r50min, r50max) &
                between(r[:,1], fmin, fmax)
            )

            if w.size > 0:
                data[ngood:ngood+w.size,:] = r[w,:]
                ngood += w.size
                nleft -= w.size

        if is_scalar:
            data=data[0,:]

        return data

    def _load_data(self):
        fname='real_galaxy_catalog_25.2_fits.fits'
        """
        fname=os.path.join(
            sys.exec_prefix,
            'share',
            'galsim',
            'COSMOS_25.2_training_sample',
            fname,
        )
        """
        fname=os.path.join(
            galsim.meta_data.share_dir,
            'COSMOS_25.2_training_sample',
            fname,
        )


        r50min,r50max=self.r50_sanity_range
        fmin,fmax=self.flux_sanity_range

        print("reading cosmos file:",fname)
        alldata=fitsio.read(fname, lower=True)
        w,=numpy.where(
            (alldata['viable_sersic']==1) &
            between(alldata['hlr'][:,0], r50min, r50max) &
            between(alldata['flux'][:,0], fmin, fmax)
        )
        print("kept %d/%d" % (w.size, alldata.size))

        self.alldata=alldata[w]

    def _make_kde(self):
        import scipy.stats

        data=numpy.zeros( (self.alldata.size, 2) )
        data[:,0] = self.alldata['hlr'][:,0]
        data[:,1] = self.alldata['flux'][:,0]

        self.kde=scipy.stats.gaussian_kde(
            data.transpose(),
            bw_method=self.kde_factor,
        )




