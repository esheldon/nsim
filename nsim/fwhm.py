from sys import stdout
import numpy

# region over which to render images and calculate likelihoods
from .sim import NSIGMA_RENDER

class Width(object):
    def __init__(self,
                 obj_model,
                 psf_model,
                 psf_T=4.0,
                 nsub=16, # for rendering
                 order=3, # for interp
                 nsub_interp=20,  # sub-interp
                 nsigma_render=NSIGMA_RENDER):

        import ngmix

        self.obj_model=obj_model
        self.psf_model=psf_model

        self.psf_T=psf_T

        self.g_prior=ngmix.priors.GPriorBA(0.3)

        self.nsub=nsub
        self.order=order

        self.nsub_interp=nsub_interp
        self.nsigma_render=nsigma_render


    def get_T_ratio_from_width_ratio(self, width_ratio, thresh, show=True,
                                     ntest=10,nrand=100):
        """
        For the input psf T and set of convolved width ratio, determine
        the required object T values

        thresh 0.5 for fwhm
        """

        minval,maxval=get_minmax_T_ratio(self.obj_model)

        T_ratios = numpy.linspace(minval, maxval, ntest)
        width_ratios = numpy.zeros(ntest)

        for i,Trat in enumerate(T_ratios):
            if (i % 5) == 0:
                stdout.write('\nTrat: %s\n' % Trat)
            else:
                stdout.write('.')

            test_T = Trat*self.psf_T

            g1,g2=self.g_prior.sample2d(nrand)
            fwhm, psf_fwhm = self._get_width_many(test_T, g1, g2, thresh)

            width_ratios[i] = (fwhm/psf_fwhm).mean()

        stdout.write('\n')

        T_ratio_a = numpy.interp(numpy.array([width_ratio]),
                                 width_ratios,
                                 T_ratios)
        T_ratio=T_ratio_a[0]


        if show:
            _plot_interp(width_ratio,
                         T_ratio,
                         width_ratios,
                         T_ratios, 
                         self.obj_model,
                         self.psf_model,
                         self.psf_T)

        return T_ratio


    def get_width(self, obj_T, obj_g1, obj_g2, thresh, show=False):
        """

        Get the full width of the convolved object and the width of the PSF

        """

        import ngmix

        obj_pars=numpy.array([0.0, 0.0,
                              obj_g1, obj_g2,
                              obj_T,
                              1.0], dtype='f8')
        psf_pars=numpy.array([0.0, 0.0,
                              0.0, 0.0,
                              self.psf_T,
                              1.0], dtype='f8')

        obj0 = ngmix.gmix.GMixModel(obj_pars, self.obj_model)
        psf  = ngmix.gmix.GMixModel(psf_pars, self.psf_model)

        obj=obj0.convolve(psf)

        T = obj.get_T()
        dims, cen=get_dims_cen(T,self.nsigma_render)

        cen += 0.1*numpy.random.randn(2)

        obj.set_cen(cen[0], cen[1])
        psf.set_cen(cen[0], cen[1])

        im = obj.make_image(dims, nsub=self.nsub)
        psf_im = psf.make_image(dims, nsub=self.nsub)

        width=measure_image_width_interp(im,
                                         thresh,
                                         nsub=self.nsub_interp,
                                         show=show,
                                         order=self.order)
        psf_width=measure_image_width_interp(psf_im,
                                             thresh,
                                             nsub=self.nsub_interp,
                                             show=show,
                                             order=self.order)

        return width, psf_width

    def _get_width_many(self, T, g1, g2, thresh):
        fwhm=numpy.zeros(g1.size)
        psf_fwhm=numpy.zeros(g1.size)
        for i in xrange(g1.size):

            fwhmi, psf_fwhmi = self.get_width(T, g1[i], g2[i], thresh)
            fwhm[i] = fwhmi
            psf_fwhm[i] = psf_fwhmi
        return fwhm, psf_fwhm

def measure_image_width_interp(image, thresh, nsub, order, show=False):
    """
    e.g. 0.5 would be the FWHM

    parameters
    ----------
    image: 2-d darray
        The image to measure
    thresh: 
        threshold is, e.g. 0.5 to get a Full Width at Half max
    nsub: int
        Number of sub-pixel
    order: order for interp
        Usually 3

    output
    ------
    width:
    """
    from numpy import array, sqrt, zeros, pi, where

    nim0 = image.copy()
    maxval=image.max()
    nim0 *= (1./maxval)

    nim = _make_interpolated_image(nim0, nsub, order, show=show)

    w=where(nim > thresh)

    area = w[0].size
    width = 2*sqrt(area/pi)/nsub

    return width


def _make_interpolated_image(im, nsub, order, show=False):
    """
    Make a new image linearly interpolated 
    on a nsubxnsub grid
    """
    # mgrid is inclusive at end when step
    # is complex and indicates number of
    # points
    import scipy.ndimage

    zimage=scipy.ndimage.zoom(im, nsub, order=order)
    if show:
        import images
        images.multiview(im)
        images.multiview(zimage)

    return zimage



def get_dims_cen(T, nsigma_render):
    """
    Based on T, get the required dimensions and a center
    """
    sigma=numpy.sqrt(T/2.)
    dims = numpy.array( [2.*sigma*nsigma_render]*2 )
    cen = numpy.array( [(dims[0]-1.)/2.]*2 )

    return dims, cen

def _plot_interp(width_ratio_best, T_ratio_best, width_ratios, T_ratios, 
                 obj_model, psf_model, psf_T):

    import biggles
    plt=biggles.FramedPlot()

    pts=biggles.Points(width_ratios, T_ratios,
                       type='filled circle',
                       color='blue')
    crv=biggles.Curve(width_ratios, T_ratios, type='solid')

    pt = biggles.Point( width_ratio_best, T_ratio_best,
                       type='filled diamond',
                       color='red',
                       size=2)

    plt.add( crv, pts, pt ) 
    plt.title='obj: %s psf: %s psf_T: %.2f' % (obj_model,psf_model,psf_T)
    plt.xlabel = r'$FWHM/FWHM_{PSF}$'
    plt.ylabel = r'$T/T_{PSF}$'

    plt.show()

def get_minmax_T_ratio(obj_model):
    if obj_model=='gauss':
        minval=0.1
        maxval=5.0
    elif obj_model=='exp':
        minval=0.2
        maxval=3.0
    elif obj_model=='dev':
        minval=2.0
        maxval=10.
    else:
        raise ValueError("unsupported model: %s" % obj_model)

    return minval, maxval
