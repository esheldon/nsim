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

    def go(self):
        import scipy.optimize
        import ngmix
        guess=self._get_guess()
        n_prior_pars=0

        res=ngmix.fitting.run_leastsq(self._errfunc,
                                      guess,
                                      n_prior_pars,
                                      maxfev=4000)

        self._result=res

    def _errfunc(self, pars):
        model = pars[0] + (self.svals - pars[1])**2/2.0/pars[2]**2

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

    res=fitter.get_result()
    return res['pars'][1], res['pars'][2]

