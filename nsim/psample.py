from __future__ import print_function
import numpy
from . import files

class PSampler(dict):
    """
    given a set of measurements and sensitvities from a deep run
    and a noisy run, calculate p(shear)
    """
    def __init__(self, run, is2n, itrial):
        self._setup(run)

        self['is2n'] = is2n
        self['itrial'] = itrial
        
        self['starti']=2
        self['nsigma']=4.5
        self._load_deep_data()
        self._set_shear_grid()

    def write(self):
        import fitsio
        fname=files.get_psample_summed_url(self['run'],
                                           self['is2n'],
                                           itrial=self['itrial'])
        print("writing:",fname)
        with fitsio.FITS(fname,'rw',clobber=True) as fits:
            fits.write(self.s1grid, extname="shear_grid")
            fits.write(self.lnp, extname="lnp")

    def go(self):
        """
        process all objects in this trial split
        """
        data=self._read_data()
        nobj=data.size

        for i in xrange(nobj):
            print("%d/%d" % (i+1,nobj))
            self._add_sums(data, i)

    def _add_sums(self, data, i):
        """
        sum over prior galaxies, evaluating at the sheared
        version of the pars, as determined by the sensitivity
        """
        from scipy.stats import multivariate_normal

        starti=self['starti']
        pars = data['pars_conv'][i,starti:]
        cov = data['pars_conv_cov'][i,starti:,starti:]

        dist = multivariate_normal(mean=pars, cov=cov)

        deep_pars=self.deep_pars
        deep_sens=self.deep_sens

        for i1,s1 in enumerate(self.s1grid):

            sheared_pars = deep_pars + deep_sens*s1
            
            likes = dist.pdf(sheared_pars)

            self.lnp[i1] += numpy.log(likes.sum()) 

    def _set_shear_grid(self):
        from .util import get_shear_grid

        self.s1grid,self.s2grid = get_shear_grid(self.sim_conf, self)

        # currently only looking at shear1
        self.lnp = numpy.zeros(self['shear_grid']['dims'][0])

    def _read_data(self):
        data=files.read_output(self['run'],
                               self['is2n'],
                               itrial=self['itrial'])
        return data


    def _load_deep_data(self):
        import fitsio
        run=self['deep_data_run']
        data=files.read_output(run, 0)

        # note g_sens (the sens from metacal) not g_sens_model
        # which means something else for max like fitting
        # make sure to convert to native data type

        si=self['starti']
        self.deep_pars=numpy.array(data['pars_noshear'][:,si:],
                                   dtype='f8',
                                   copy=True)
        self.deep_sens=numpy.array(data['sens'][:,si:],
                                   dtype='f8',
                                   copy=True)


    def _setup(self, run, **keys):
        """
        Check and set the configurations
        """

        run_conf = files.read_config(run)
        sim_conf = files.read_config(run_conf['sim'])

        if sim_conf['name'] != run_conf['sim']:
            err="sim name in run config '%s' doesn't match sim name '%s'"
            raise ValueError(err % (run_conf['sim'],sim_conf['name']))

        self.sim_conf=sim_conf
        self.update(run_conf)
