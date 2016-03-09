import numpy

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
