from __future__ import print_function
import numpy

def ring_select(logic):

    weven = numpy.arange(0,logic.size,2)
    wodd  = numpy.arange(1,logic.size,2)

    wboth, = numpy.where(  logic[wodd] & logic[weven] )

    w=numpy.concatenate( (wodd[wboth], weven[wboth] ) )

    w.sort()
    return w

def get_shear_grid(simconf, runconf):
    gc = runconf['shear_grid']

    desired_err=runconf['desired_err']
    true_shear = simconf['shear'][0]
    nsigma = gc['nsigma']

    shmin = true_shear - nsigma*desired_err
    shmax = true_shear + nsigma*desired_err
    shear_grid = numpy.linspace(shmin,
                                shmax,
                                gc['npoints'])
    return shear_grid

class lnp_fitter(object):
    def __init__(self, svals, lnp):
        self.svals = svals
        self.lnp = lnp
        self.mlnp = -lnp

        self._set_guesses()

    def get_result(self):
        return self._result

    def go(self, guess=None):
        import scipy.optimize
        import ngmix

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


def fit_lnp_shear(simconf, runconf, lnp_shear):
    """
    fit -lnp_shear to a parabola
       model = a + ( (shear - shear_mean)^2/2.0/shear_err^2 )
    free parameters are a, shear_mean, shear_err
    """
    from nsim.util import get_shear_grid

    s=get_shear_grid(simconf, runconf)

    fitter=lnp_fitter(s, lnp_shear)
    fitter.go()

    return fitter

