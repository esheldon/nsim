from __future__ import print_function

try:
    xrange
except:
    xrange=range

import os, sys
import logging
import numpy
from numpy import exp, log, zeros, ones, sqrt, newaxis
import ngmix

import esutil as eu
from esutil.random import srandu
from ngmix.gexceptions import GMixRangeError

from . import files

logger = logging.getLogger(__name__)

class TryAgainError(Exception):
    """
    signal to skip this image(s) and try a new one
    """
    def __init__(self, message):

        # Call the base class constructor with the parameters it needs
        Exception.__init__(self, message)

def setup_logging(level):
    if level=='info':
        l=logging.INFO
    elif level=='debug':
        l=logging.DEBUG
    elif level=='warning':
        l=logging.WARNING
    elif level=='error':
        l=logging.ERROR
    else:
        l=logging.CRITICAL

    logging.basicConfig(stream=sys.stdout, level=l)

def log_pars(pars, fmt='%8.3g',front=None):
    """
    print the parameters with a uniform width
    """

    s = []
    if front is not None:
        s.append(front)
    if pars is not None:
        fmt = ' '.join( [fmt+' ']*len(pars) )
        s.append( fmt % tuple(pars) )
    s = ' '.join(s)

    logger.debug(s)

class Namer(object):
    """
    create strings with a specified front prefix
    """
    def __init__(self, front=None, back=None):
        if front=='':
            front=None
        if back=='' or back=='noshear':
            back=None

        self.front=front
        self.back=back

        if self.front is None and self.back is None:
            self.nomod=True
        else:
            self.nomod=False



    def __call__(self, name):
        n = name
        if not self.nomod:
            if self.front is not None:
                n = '%s_%s' % (self.front, n)
            if self.back is not None:
                n = '%s_%s' % (n, self.back)
        
        return n

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

def get_true_shear(conf, index=0):
    """
    if shear is constant, return that, otherwise sample it
    """

    shc=conf['shear']
    if isinstance(shc, dict):
        if 'shears' in shc:
            shear=shc['shears'][index]
        else:
            shear=shc['mean']
    else:
        shear=shc

    return shear


def write_fits(filename, data, header=None):
    """
    Assume condor where cwd is scratch dir

    Write to cwd assuming scratch the move the file.

    The move may fail; it is retried a few times.
    """
    import fitsio
    import time

    output_file=os.path.abspath(filename)

    local_file=os.path.abspath( os.path.basename(output_file) )
    logger.info("writing local file: %s" % local_file)

    fitsio.write(local_file, data, header=header, clobber=True)

    if local_file==output_file:
        return True

    # remove if it exists
    try:
        os.remove(output_file)
    except:
        pass

    # try a few times
    logger.info("moving to: %s" % output_file)
    cmd='mv %s %s' % (local_file, output_file)
    for i in xrange(5):
        stat=os.system(cmd)
        if stat==0:
            logger.info('success')
            success=True
            break
        else:
            logger.info('error moving file, trying again')
            time.sleep(5)
            success=False

    return success

def get_weights(data,SN=0.24,field='g_cov',type='noise'):
    """
    Using Ts2n weights seems to  help when there is a mixture
    of exp and dev, essentially downweighting dev
    """

    if type=='noise':
        #print('field:',field,'SN:',SN)
        denom = (
            2*SN**2
            + data[field][:,0,0]
            + 2*data[field][:,0,1]
            + data[field][:,1,1]
        )
        wts = 1.0/denom

    elif type == 's2n':
        # fit to run-bd06zmcal-degrade03 
        fitcoeff=numpy.array([ 0.05676219, -0.36838486,  0.92523771, -1.07835695,  0.49725526])
        ply=numpy.poly1d(fitcoeff)

        #print('field:',field,'SN:',SN)
        #print(ply)

        logs2n = numpy.log10(data[field])
        var_per_component = ply(logs2n)

        #denom = ( 2*SN**2 + 2*var_per_component )
        denom = ( 2*SN**2 + 2*var_per_component*2 )
        wts = 1.0/denom

    elif type=='s2n_simple':

        err_per_component = 1.0/data[field]
        denom = ( 2*SN**2 + 2*err_per_component**2 )
        wts = 1.0/denom

    elif type is None:
        wts = ones(data.size)

    else:
        raise ValueError("bad weight type: '%s'" % type)

    wts *= (1.0/wts.max())
    wfrac=wts.sum()/wts.size
    print("    weighted frac: %.3g" % wfrac)
    return wts


def fit_prior(run, is2n=0, field='pars_noshear',show=False):
    import biggles
    import esutil as eu
    from sklearn.mixture import GMM
    import fitsio
    gm=ngmix.gmix.GMixND()

    alldata=files.read_output(run, is2n)
    data=alldata[field]

    # fit |g|
    g=sqrt(data[:,2]**2 + data[:,3]**2)

    gp=ngmix.priors.GPriorBA()

    bs=eu.stat.Binner(g)
    bs.dohist(nbin=100)
    bs.calc_stats()
    xvals=bs['center']
    yvals=bs['hist'].astype('f8')
    gp.dofit(xvals, yvals)
    rg=gp.sample1d(1000000)

    #eta1,eta2,good=ngmix.shape.g1g2_to_eta1eta2_array(data[:,2],
    #                                                  data[:,3])
    #eta = 2*numpy.arctanh(g)
    #logg=log(g)


    outfile=files.get_fitprior_url(run, is2n)
    epsfile=files.get_fitprior_url(run, is2n, ext='eps')

    # fit TF with n-dimensional gaussian
    ngauss=20
    gm.fit(data[:, 4:], ngauss, min_covar=1.0e-4)

    '''
    print("fitting log(g)")
    g_gmm=GMM(n_components=10,
                 n_iter=5000,
                 min_covar=1.0e-4,
                 covariance_type='full')
    g_gmm.fit(eta[:,newaxis])

    if not g_gmm.converged_:
        print("DID NOT CONVERGE")

    rg=g_gmm.sample(1000000)
    gplt=biggles.plot_hist(eta,nbin=100,xlabel='|g|',
                           norm=1,
                           visible=False)
    '''
    gplt=biggles.plot_hist(g,nbin=100,xlabel='|g|',
                           norm=1,
                           visible=False)
    #gplt=biggles.plot_hist(rg[:,0],nbin=100,
    gplt=biggles.plot_hist(rg,nbin=100,
                           plt=gplt,
                           norm=1,
                           color='red',
                           visible=False)

    print("saving mixture to:",outfile)
    gm.save_mixture(outfile)
    with fitsio.FITS(outfile,'rw') as fits:
        gout=numpy.zeros(1, dtype=[('sigma','f8')])
        gout['sigma'] = gp.fit_pars[1]

        fits.write(gout, extname='gfit')

    #
    # plots
    #
    r=gm.sample(1000000)

    tab=biggles.Table(2,2)

    Tplt=biggles.plot_hist(data[:,4], nbin=100,
                          xlabel='log(T)',
                          norm=1,
                          visible=False)
    biggles.plot_hist(r[:,0], nbin=100,
                      color='red',
                      norm=1,
                      plt=Tplt,
                      visible=False)

    Fplt=biggles.plot_hist(data[:,5], nbin=100,
                          xlabel='log(F)',
                          norm=1,
                          visible=False)
    biggles.plot_hist(r[:,1], nbin=100,
                      color='red',
                      norm=1,
                      plt=Fplt,
                      visible=False)

    '''
    eta1plt=biggles.plot_hist(eta1, nbin=100,
                              xlabel='eta1',
                              norm=1,
                              visible=False)
    eta2plt=biggles.plot_hist(eta2, nbin=100,
                              xlabel='eta2',
                              norm=1,
                              visible=False)
    '''

    tab[0,0] = Tplt
    tab[0,1] = Fplt
    #tab[1,0] = eta1plt
    #tab[1,1] = eta2plt
    tab[1,0] = gplt
    tab[1,1] = gplt

    tab.aspect_ratio=0.5

    print("writing:",epsfile)
    tab.write_eps(epsfile)

    if show:
        biggles.configure('screen','width',1200)
        biggles.configure('screen','height',700)
        tab.show()

def load_gmixnd(spec, rng=None):
    if 'means' in spec:
        pdf=ngmix.gmix.GMixND(
            weights=spec['weights'],
            means=spec['means'],
            covars=spec['covars'],
            rng=rng,
        )
    else:
        if 'run' in spec:
            extra=spec['extra']
            fname=files.get_fitprior_url(spec['run'], 0, extra=extra)
            
        else:
            fname=files.get_extra_url(spec['file'])

        pdf=ngmix.gmix.GMixND(file=fname, rng=rng)

    if 'cov_factor' in spec:
        logger.info("    using cov factor: %s" % spec['cov_factor'])
        pdf.covars *= spec['cov_factor']

    return pdf


def get_sheared_jacobian(jacin, shear):

    jac=jacin.copy()

    fac = 1.0/sqrt(1.-shear.g1**2 - shear.g2**2)

    S = numpy.array([[ 1.+shear.g1,      shear.g2  ],
                     [    shear.g2  , 1.-shear.g1  ]])
    S *= fac

    dat=jac._data
    jmat = numpy.array([[dat['dudrow'][0], dat['dudcol']],
                        [dat['dvdrow'][0], dat['dvdcol']]])

    new_jmat = numpy.dot(S, jmat)


    dat['dudrow'][0] = new_jmat[0,0]
    dat['dudcol'][0] = new_jmat[0,1]
    dat['dvdrow'][0] = new_jmat[1,0]
    dat['dvdcol'][0] = new_jmat[1,1]

    return jac

def find_kmax(kr, ki, min_rel_val):
    """
    get the maximum radius in k space for which the value is larger
    than the indicated value relative to the maximum

    """
    import deconv

    ps = kr.array**2 + ki.array**2
    amp = numpy.sqrt(ps)

    dims=amp.shape
    if (dims[0] % 2) == 0:
        rowshift=0.5
    else:
        rowshift=0.0

    if (dims[1] % 2) == 0:
        colshift=0.5
    else:
        colshift=0.0

    cen=(numpy.array(dims)-1.0)/2.0 + [rowshift,colshift]

    rows,cols=deconv.util.make_rows_cols(
        dims,
        cen=cen,
    )

    r2 = rows**2 + cols**2

    maxval = amp.max()
    minval = maxval*min_rel_val**2

    w=numpy.where(amp > minval)
    if w[0].size == 0:
        raise DeconvRangeError("no good psf values in k space")

    kmax = numpy.sqrt(r2[w].max())
    return kmax

def get_shrink_factor(shear):
    e1,e2=ngmix.shape.g1g2_to_e1e2(shear.g1,shear.g2)
    e=sqrt(e1**2 + e2**2)

    if e >= 0.99999:
        raise GMixRangeError("e out of bounds")
    s = ( (1+e)/(1-e) )**0.25
    return s

def complex_multiply(a, b, c, d, scratch, real_res, imag_res):
    """
    (a + i *b) * (c + i *d)
    =
    (ac-bd) + i* (ad + bc)
    """

    numpy.multiply(a, c, real_res)
    numpy.multiply(b, d, scratch)

    real_res -= scratch

    numpy.multiply(a, d, imag_res)
    numpy.multiply(b, c, scratch)

    imag_res += scratch

def get_kmom_shear(data, step=0.01):
    h,rev=eu.stat.histogram(data['shear_index'], rev=True)

    nbin=h.size

    dt=[
        ('shear_true','f8',2),
        ('shear','f8',2),
        ('shear_err','f8',2),
    ]
    out=numpy.zeros(nbin, dtype=dt)
    for i in xrange(nbin):
        sys.stdout.write("."); sys.stdout.flush()
        w=rev[ rev[i]:rev[i+1] ]
        R1, R2 = get_kmom_R(data[w], step=step)
        R1mean=R1.mean()
        R2mean=R2.mean()

        R=numpy.array([R1mean, R2mean])

        g=data['mcal_g'][w].mean(axis=0)
        gerr=data['mcal_g'][w].std(axis=0)/sqrt(w.size)

        out['shear'][i,:] = g/R
        out['shear_err'][i,:] = gerr/R

        out['shear_true'][i,:] = data['shear_true'][w[0]]

    print()
    return out


def get_kmom_R(data, step=0.01):
    """
    need selections
    """

    e1 = data['mcal_g'][:, 0]
    e2 = data['mcal_g'][:, 1]

    # not normalized
    T  = data['mcal_pars'][:, 4]#/data['mcal_wsum']

    M1R, TR1, M2R, TR2 = _get_kmom_R(data, 'mcal_pars_', step=step)
    wM1R, wTR1, wM2R, wTR2 = _get_kmom_R(data, 'mcal_pars_w', step=step)

    R1tot = M1R/T - e1/T*TR1
    R2tot = M2R/T - e2/T*TR2

    wR1 = wM1R/T - e1/T*wTR1
    wR2 = wM2R/T - e2/T*wTR2

    R1 = R1tot - wR1
    R2 = R2tot - wR2

    return R1, R2

def _get_kmom_R(data, front, step=0.01):
    pars_1p=data['%s%s' % (front,'1p')]
    pars_1m=data['%s%s' % (front,'1m')]
    pars_2p=data['%s%s' % (front,'2p')]
    pars_2m=data['%s%s' % (front,'2m')]

    M1R = (pars_1p[:,2]-pars_1m[:,2])/(2*step)
    TR1 = (pars_1p[:,4]-pars_1m[:,4])/(2*step)

    M2R = (pars_2p[:,3]-pars_2m[:,3])/(2*step)
    TR2 = (pars_2p[:,4]-pars_2m[:,4])/(2*step)

    return M1R, TR1, M2R, TR2


def fit_line(y, x):
    """
    fit a line such that y = ax + b
        
    We can rewrite the line equation as

        gamma = Ap 

    where A = [[x 1]], and the ones are inserted for the intercept

    this is accomplished using vstack

    A = vstack( [x, ones(len(x))] ).T

    coeff,_,_,_ = np.linalg.lstsq(A, y)[0]

    a,b = coeff
    """
    from numpy import vstack, ones

    A = vstack([x, ones(len(x))]).T

    coeff, _, _, _ = numpy.linalg.lstsq(A, y)
    a,b = coeff
    return a,b


def convert_run_to_seed(run):
    """
    convert the input config file name to an integer for use
    as a seed
    """
    import hashlib

    h = hashlib.sha256(run.encode('utf-8')).hexdigest()
    seed = int(h, base=16) % 2**30 

    logger.info("got seed %d from run %s" % (seed,run))

    return seed

