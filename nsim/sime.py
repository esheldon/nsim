from __future__ import print_function

from pprint import pprint
import numpy

from . import files
from . import psfs
from . import objects
from . import observations

from .shearpdf import get_shear_pdf


class Sim(dict):
    def __init__(self, sim_conf, seed):
        import galsim

        self._load_config(sim_conf)

        print("using seed:",seed)

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


    def _load_config(self, confin):

        if isinstance(confin,dict):
            # full config dictionary was input
            conf=confin
        else:
            if '.yaml' in confin:
                # full path given
                conf=files.read_yaml(confin)
            else:
                # identifier given, assumed to be
                # in the "usual" place
                conf=files.read_config(confin)
            # in this case, we offer to set the seed
            # if it is not in the file

        self.update(conf)
