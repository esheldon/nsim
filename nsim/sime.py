from __future__ import print_function
try:
    xrange
except:
    xrange=range
try:
    raw_input
except:
    raw_input=input


import os, sys
from pprint import pprint

import numpy
import fitsio

import ngmix
import esutil as eu
from esutil.numpy_util import between

from . import psfs
from . import objects
from . import sim as ngmixsim
from ngmix.priors import srandu

from .util import TryAgainError, load_gmixnd

from . import pdfs
from .pdfs import DiscreteSampler, PowerLaw
import galsim


class SimE(dict):
    def __init__(self, sim_conf):
        self.update(sim_conf)

        seed=self.get('seed',None)
        print("    using seed:",self['seed'])

        # seeding both the global and the local rng.  With the
        # local, we produce the same sim independent of the fitting
        # code which may use the global.
        numpy.random.seed(seed)
        self.rng=numpy.random.RandomState(seed=numpy.random.randint(0,2**30))
        self.galsim_rng = galsim.BaseDeviate(numpy.random.randint(0,2**30))

        self._setup()

        pprint(self)

        self._set_makers()

    def _set_makers(self):
        psf_maker = psfs.get_psf_maker(self['psf'],self.rng)
        self._object_maker = objects.get_object_maker(self['obj_model'], psf_maker,self.rng)
