import numpy
from ngmix import Shape
from esutil.stat import print_stats

def get_shear_pdf(conf, rng=None):
    from .shearpdf import ConstShearSelector, ConstShearGenerator

    if 'shear' in conf:
        shconf = conf['shear']
        # shears are imbedded in the config
        if shconf['type'] == 'const':
            pdf = ConstShearSelector(shconf['shears'], rng=rng)
        elif shconf['type'] == 'const-dist':
            # a seed is specified and we generate them
            rng=numpy.random.RandomState(seed=shconf['seed'])
            pdf = ConstShearGenerator(
                rng,
                shconf['nshear'],
                min_shear=shconf['min_shear'],
                max_shear=shconf['max_shear'],
            )
        else:
            raise ValueError("only shear 'const' for now")

    else:
        pdf=None

    return pdf

class ShearGeneratorBase(object):
    def get_shear(self):
        raise NotImplementedError("implement get_shear()")

class ConstShearSelector(ShearGeneratorBase):
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

class ConstShearGenerator(ShearGeneratorBase):
    def __init__(self, rng, nshear, min_shear=0.1, max_shear=0.08):

        self.rng=rng
        self.nshear=nshear
        self.min_shear=min_shear
        self.max_shear=max_shear

        self.gen_shear()

    def gen_shear(self):
        g = self.rng.uniform(
            low=self.min_shear,
            high=self.max_shear,
            size=self.nshear,
        )
        theta = self.rng.uniform(
            low=0.0,
            high=numpy.pi*2,
            size=self.nshear,
        )

        g1=g*numpy.cos(2.0*theta)
        g2=g*numpy.sin(2.0*theta)

        print("generated shears:")
        print_stats(g1)
        print_stats(g2)

        shears=[]
        for i in xrange(self.nshear):
            shears.append( Shape(g1[i], g2[i]) )

        self.shears=shears

    def get_shear(self):
        """
        return a random shear from the input list, plus an index
        """
        ri = self.rng.randint(0, self.nshear)
        return self.shears[ri], ri


