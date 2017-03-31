from __future__ import print_function

import os, sys
from pprint import pprint

import numpy
import fitsio

import ngmix
import esutil as eu
from esutil.numpy_util import between

from . import psfs
from . import objects
from . import observations

from .util import TryAgainError

from .shearpdf import get_shear_pdf

def get_sim(sim_conf):
    return Sim(sim_conf)

class Sim(dict):
    def __init__(self, sim_conf):
        import galsim
        self.update(sim_conf)

        seed=self['seed']
        print("using seed:",self['seed'])

        # seeding both the global and the local rng.  With the
        # local, we produce the same sim independent of the fitting
        # code which may use the global.
        numpy.random.seed(seed)
        self.rng=numpy.random.RandomState(seed=numpy.random.randint(0,2**30))
        self.galsim_rng = galsim.BaseDeviate(self.rng.randint(0,2**30))

        pprint(self)

        self._set_makers()

    def __call__(self):
        return self._image_maker()

    def _set_makers(self):
        psf_maker    = psfs.get_psf_maker(self['psf'], self.rng)

        # for multi-band, we will make multiple of these
        object_maker = objects.get_object_maker(
            self['object'],
            self.rng,
            self.galsim_rng,
        )

        if 'shear' in self:
            shear_pdf = get_shear_pdf(self['shear'], self.rng)
        else:
            shear_pdf = None

        self._image_maker  = observations.get_observation_maker(
            self['images'],
            psf_maker,
            object_maker,
            self.rng,
            self.galsim_rng,
            shear_pdf=shear_pdf,
        )


