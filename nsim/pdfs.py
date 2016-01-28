import numpy

class DiscreteSampler(object):
    def __init__(self, vals):
        self.vals=vals
        self.nvals=len(vals)

    def sample(self):
        ri=numpy.random.randint(0, self.nvals)
        return self.vals[ri]
