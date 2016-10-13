import numpy
import esutil as eu

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

        hlr=self.data['hlr'][ri,0]
        flux=self.data['flux'][ri,0]

        return hlr, flux


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
