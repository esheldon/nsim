from __future__ import print_function
import numpy
from numpy import exp, zeros, sqrt
import ngmix
from esutil.random import srandu

def ring_select(logic):

    weven = numpy.arange(0,logic.size,2)
    wodd  = numpy.arange(1,logic.size,2)

    wboth, = numpy.where(  logic[wodd] & logic[weven] )

    w=numpy.concatenate( (wodd[wboth], weven[wboth] ) )

    w.sort()
    return w

def get_shear_grid(simconf, runconf, get_mesh=False):
    """
    get s1 and s2 grids

    parameters
    ----------
    simconf: dict
        config for sim, from which we get the true shear
    runconf: dict
        run config, from which we get grid parameters
    get_mesh: bool
        If True, return a mesh grid rather than 1-d grids
        in each dimension
    """
    gc = runconf['shear_grid']

    desired_err=runconf['desired_err']
    if isinstance(simconf['shear'], dict):
        true_shear=simconf['shear']['mean']
    else:
        true_shear=simconf['shear']

    nsigma = gc['nsigma']

    sh1min = true_shear[0] - nsigma*desired_err
    sh1max = true_shear[0] + nsigma*desired_err
    sh2min = true_shear[1] - nsigma*desired_err
    sh2max = true_shear[1] + nsigma*desired_err

    g1g = numpy.linspace(sh1min,
                         sh1max,
                         gc['dims'][0])
    g2g = numpy.linspace(sh2min,
                         sh2max,
                         gc['dims'][1])

    if get_mesh:
        grid1, grid2 = numpy.meshgrid(g1g, g2g)
        return grid1, grid2
    else:
        return g1g, g2g

class lnp_fitter2d(object):
    def __init__(self, s1vals, s2vals, lnp):
        self.s1vals = s1vals
        self.s2vals = s2vals
        self.lnp = lnp
        self.p=exp(lnp-lnp.max())

        self.scale1=s1vals[1]-s1vals[0]
        self.scale2=s2vals[1]-s2vals[0]

        #print("scale: (%g,%g)" % (self.scale1,self.scale2))
        self.jacobian=ngmix.Jacobian(0.0,
                                     0.0,
                                     self.scale1,
                                     0.0,
                                     0.0,
                                     self.scale2)


        self._set_guesses()

    def get_result(self):
        return self._result

    def go(self, guess=None):

        if guess is None:
            guess=self._get_guess()
        else:
            guess=self._transform_guess(guess)

        ngmix.print_pars(guess,front="    guess: ")

        obs=ngmix.Observation(self.p, jacobian=self.jacobian)

        fitter=ngmix.fitting.LMSimple(obs, 'gauss')
        fitter.go(guess)

        res=fitter.get_result()

        ngmix.print_pars(res['pars'],front="    pars:  ")

        self._add_shear_info(res)


        self._result=res
        self._fitter=fitter

    def _transform_guess(self, guess):
        # s1,s2, se1,se2, T, Amp
        nguess=zeros(6)

        Irr = guess['shear_cov'][0,0]
        Irc = 0.0
        Icc = guess['shear_cov'][1,1]

        T = Irr + Icc

        Flux = self.p.sum()*self.scale1*self.scale2

        # guess the middle
        nguess[0] = (self.s1vals[-1]-self.s1vals[0])/2.
        nguess[1] = (self.s2vals[-1]-self.s2vals[0])/2.
        nguess[2] = 0.0 + 0.01*srandu()
        nguess[3] = 0.0 + 0.01*srandu()
        nguess[4] = T
        nguess[5] = Flux

        return nguess
        

    def _add_shear_info(self, res):
        """
        the shear is the centroid

        the shear_cov is the [Ixx,Ixy,
                              Ixy,Iyy]
        """

        pars=res['pars']

        # these are not the shear, but to get the covariance
        # matrix of the gaussian, which incodes the covarianc
        # matrix of the shear
        g1,g2=pars[2:2+2]
        T=pars[4]

        e1,e2=ngmix.shape.g1g2_to_e1e2(g1,g2)

        Irr = (T/2)*(1-e1)
        Irc = (T/2)*e2
        Icc = (T/2)*(1+e1)

        res['shear'] = pars[0:0+2].copy()
        res['shear'][0] += self.s1vals[0]
        res['shear'][1] += self.s2vals[0]

        res['shear_cov'] = numpy.array( [ [Irr, Irc],
                                          [Irc, Icc] ], dtype='f8' )

    def doplot(self, show=False):
        import images

        gm=self._fitter.get_gmix()

        model=gm.make_image(self.p.shape, jacobian=self.jacobian)

        plt=images.compare_images(self.p, model,
                                  label1='p(g1,g2)',
                                  label2='model',
                                  show=show)

        return plt

    def _set_guesses(self):
        self.a_guess = self.p.max()
        self.sh_guess = [self.s1vals.mean(), self.s2vals.mean()]

        # we ususally have use a range of ~+/-4 sigma around the mean
        # so sigma should be [smax-smin]/2./4.0
        sigma_guess = 0.5*(self.s1vals.std() + self.s2vals.std())/2./4.
        self.Tguess= 2*sigma_guess**2

    def _get_guess(self):
        return numpy.array([self.sh_guess[0] + 0.01*srandu(),
                            self.sh_guess[1] + 0.01*srandu(),
                            0.0+0.01*srandu(),
                            0.0+0.01*srandu(),
                            self.Tguess*(1.0+0.01*srandu()),
                            self.a_guess*(1.0+0.01*srandu())
                            ])



class lnp_fitter1d(object):
    def __init__(self, svals, lnp):
        self.svals = svals
        self.lnp = lnp
        self.mlnp = -lnp

        self._set_guesses()

    def get_result(self):
        return self._result

    def go(self, guess=None):
        import scipy.optimize

        if guess is None:
            guess=self._get_guess()

        n_prior_pars=0

        res=ngmix.fitting.run_leastsq(self._errfunc,
                                      guess,
                                      n_prior_pars,
                                      maxfev=4000)

        self._result=res

    def doplot(self, show=False):
        import biggles

        plt=biggles.FramedPlot()
        plt.aspect_ratio=1
        plt.xlabel='shear'
        plt.ylabel='-ln(prob)'

        pts=biggles.Points(self.svals, self.mlnp, type='filled circle')

        pars=self._result['pars']
        model=self.get_model(pars)

        c=biggles.Curve(self.svals, model, color='blue')

        pts.label='data'
        c.label='model'

        key=biggles.PlotKey(0.9, 0.9, [pts,c],
                            halign='right')

        plt.add(pts, c, key)

        if show:
            plt.show()
        return plt

    def get_model(self, pars):
        model = pars[0] + (self.svals - pars[1])**2/2.0/pars[2]**2
        return model

    def _errfunc(self, pars):
        model = self.get_model(pars)
        return model-self.mlnp

    def _set_guesses(self):
        self.a_guess = self.mlnp.min()
        self.sh_guess = self.svals.mean()
        # we ususally have use a range of ~+/-4 sigma around the mean
        # so sigma should be [smax-smin]/2./4.0
        self.sherr_guess = self.svals.std()/2./4.

    def _get_guess(self):
        return numpy.array([self.a_guess,
                      self.sh_guess,
                      self.sherr_guess], dtype='f8')


def fit_lnp_shear2d(simconf, runconf, lnp_shear, guess=None):
    """
    fit p to 2d gaussian

    guess is on cov {'shear_cov':[[s11,s12,
                                   s12,s22]]}

    """

    s1grid,s2grid=get_shear_grid(simconf, runconf)

    fitter=lnp_fitter2d(s1grid, s2grid, lnp_shear)
    fitter.go(guess=guess)

    return fitter

def fit_lnp_shear1d(simconf, runconf, lnp_shear):
    """
    fit -lnp_shear to a parabola
       model = a + ( (shear - shear_mean)^2/2.0/shear_err^2 )
    free parameters are a, shear_mean, shear_err
    """
    from nsim.util import get_shear_grid

    s=get_shear_grid(simconf, runconf)

    fitter=lnp_fitter1d(s, lnp_shear)
    fitter.go()

    return fitter

def get_true_shear(conf):
    """
    if shear is constant, return that, otherwise sample it
    """
    if isinstance(conf['shear'], dict):
        shear=conf['shear']['mean']
    else:
        shear=conf['shear']

    return shear
