import numpy
from ngmix import Shape

class ShearGeneratorBase(object):
    def get_shear(self):
        raise NotImplementedError("implement get_shear()")

class ConstShearGenerator(ShearGeneratorBase):
    def __init__(self, shears, rng=None):

        if rng is None:
            rng=numpy.random.RandomState()
        self.rng=rng

        if not isinstance(shears[0], list):
            shears = [shears]

        shears = [Shape(s[0],s[1]) for s in shears]

        self.nshear=len(shears)
        self.shears=shears

    def get_shear(self):
        """
        return a random shear from the input list, plus an index
        """
        ri = self.rng.randint(0, self.nshear)
        return self.shears[ri], ri
