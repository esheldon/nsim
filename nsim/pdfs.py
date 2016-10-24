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
        fname=os.path.join(
            sys.exec_prefix,
            'share',
            'galsim',
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



class RandomWalkGalaxy(object):
    """
    place point sources randomly by letting them do
    a random walk
    
    paramters
    ---------
    hlr: float
        Half-light radius for the final distribution
    flux: float
        Total flux of the set of points
    npoints: int, optional
        Total number of points to use.  Default 100
    nstep: int, optional
        Number of steps in random walk.  Default 40
    """
    def __init__(self, hlr, flux, npoints=100, nstep=40):

        self.hlr=hlr
        self.flux=flux
        self.npoints=npoints
        self.nstep=nstep

        self.factor = numpy.sqrt(nstep)/2.09
        self.scale = hlr/self.factor

        self._set_points()
        self._set_gsobj()

    def get_gsobj(self):
        """
        get the galsim object
        """
        return self.g

    def _set_gsobj(self):

        fluxper=self.flux/self.npoints
        gaussians=[]

        points=self.points
        for i in xrange(points.shape[0]):
            dx,dy = points[i]

            g=galsim.Gaussian(sigma=1.0e-3, flux=fluxper)
            g = g.shift(dx=dx, dy=dy)

            gaussians.append(g)

        self.g = galsim.Add(gaussians)

    def _set_points(self):

        scale=self.scale
        npoints=self.npoints
        nstep=self.nstep

        pts=numpy.zeros( (npoints, 2) )

        for i in xrange(npoints):
            x=0.0
            y=0.0

            for istep in xrange(nstep):

                r = scale*numpy.random.random()
                angle = 2*numpy.pi*numpy.random.random()

                dx = r*numpy.cos(angle)
                dy = r*numpy.sin(angle)

                x += dx
                y += dy

            pts[i,0] = x
            pts[i,1] = y

        self.points=pts

