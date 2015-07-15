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

            if (i % 100) == 0:
                self._plot_shear()

            self._add_sums(data, i)

    def _add_sums(self, data, i):
        """
        sum over prior galaxies, evaluating at the sheared
        version of the pars, as determined by the sensitivity
        """
        from scipy.stats import multivariate_normal

        #starti=self['starti']
        #fname,cname='pars_conv','pars_conv_cov'
        fname,cname='pars','pcov'
        pars = data[fname][i]
        cov = data[cname][i]

        dist = multivariate_normal(mean=pars, cov=cov)

        deep_pars=self.deep_pars
        deep_sens=self.deep_sens

        sheared_pars=deep_pars.copy()

        if True:
            self._compare_dist(dist)

        for i1,s1 in enumerate(self.s1grid):

            #sheared_pars[:] = deep_pars
            sheared_pars[:,2] = deep_pars[:,2]
            sheared_pars[:,2] += deep_sens[:,2]*s1
            
            likes = dist.pdf(sheared_pars)

            likesum = likes.sum()

            effnum=likesum/likes.max()
            print("effnum:",effnum)
            if likesum > 0:
                self.lnp[i1] += numpy.log( likesum )

    def _set_shear_grid(self):
        from .util import get_shear_grid

        shconf=self['shear_grid']
        if 's1range' in shconf:
            rng1=shconf['s1range']
            rng2=shconf['s2range']

            self.s1grid=numpy.linspace(rng1[0],
                                       rng1[1],
                                       shconf['dims'][0])
            self.s2grid=numpy.linspace(rng2[0],
                                       rng2[1],
                                       shconf['dims'][1])
        else:
            self.s1grid,self.s2grid = get_shear_grid(self.sim_conf, self)

        # currently only looking at shear1
        self.lnp = numpy.zeros(self['shear_grid']['dims'][0])

    def _plot_shear(self, show=False):
        import biggles
        from biggles import plot

        tab=biggles.Table(2,1)

        tab[0,0]=plot(self.s1grid,
                      self.lnp,
                      xlabel='g1',
                      ylabel='log(p)',
                      visible=False)
        p=numpy.exp(self.lnp-self.lnp.max())
        tab[1,0]=plot(self.s1grid,
                      p,
                      xlabel='g1',
                      ylabel='p',
                      visible=False)

        pngfile=files.get_psample_summed_url(self['run'],
                                             self['is2n'], 
                                             itrial=self['itrial'],
                                             ext='png')
        print("writing:",pngfile)
        tab.write_img(800,800,pngfile)

        if show:
            tab.show()
        return tab

    def _compare_dist(self, dist):
        import biggles
        import esutil as eu
        from biggles import plot_hist

        r=dist.rvs(self.deep_pars.shape[0])

        names=['c1','c2','M1','M2','T','I']
        grid=eu.plotting.Grid(6)
        tab=biggles.Table(grid.nrow, grid.ncol)
        for i in xrange(6):
            row,col=grid(i)

            plt=plot_hist(r[:,i], nbin=100, 
                          xlabel=names[i],
                          visible=False)
            plt=plot_hist(self.deep_pars[:,i], nbin=100, color='red', plt=plt,
                         visible=False)

            tab[row,col] = plt

        pngfile=files.get_psample_summed_url(self['run'],
                                             self['is2n'], 
                                             itrial=self['itrial'],
                                             ext='png')
 
        pngfile=pngfile.replace('.png','-hist.png')
        print(pngfile)
        tab.write_img(1200,800,pngfile)
        key=raw_input('hit a key: ')
        if key=='q':
            stop

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

        fname='pars'
        self.deep_pars=numpy.array(data[fname],
                                   dtype='f8',
                                   copy=True)
        self.deep_sens=numpy.array(data['sens'],
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
