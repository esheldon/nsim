import numpy

# region over which to render images and calculate likelihoods
from .sim import NSIGMA_RENDER

def get_width(obj_model,
              obj_T,
              psf_model,
              psf_T,
              thresh,
              nsigma_render=NSIGMA_RENDER,
              nsub=16,
              smooth=0.1):
    """

    Get the full width of the convolved object and the width of the PSF

    """

    import ngmix
    

    obj_pars=numpy.array([0.0, 0.0,
                          0.0, 0.0,
                          obj_T,
                          1.0], dtype='f8')
    psf_pars=numpy.array([0.0, 0.0,
                          0.0, 0.0,
                          psf_T,
                          1.0], dtype='f8')

    obj0 = ngmix.gmix.GMixModel(obj_pars, obj_model)
    psf  = ngmix.gmix.GMixModel(psf_pars, psf_model)

    obj=obj0.convolve(psf)

    T = obj.get_T()
    dims, cen=get_dims_cen(T,nsigma_render)

    obj.set_cen(cen[0], cen[1])
    psf.set_cen(cen[0], cen[1])

    im = obj.make_image(dims, nsub=nsub)
    psf_im = psf.make_image(dims, nsub=nsub)

    width=measure_image_width_erf(im, thresh, smooth)
    psf_width=measure_image_width_erf(psf_im, thresh, smooth)

    return width, psf_width


def measure_image_width_erf(image, thresh, smooth):
    """
    Measure width at the given threshold using an erf to smooth the contour.
    
    e.g. 0.5 would be the FWHM

    parameters
    ----------
    image: 2-d darray
        The image to measure
    thresh: number
        threshold is, e.g. 0.5 to get a Full Width at Half max
    smooth: float
        The smoothing scale for the erf.  This should be between 0 and 1. If
        you have noisy data, you might set this to the noise value or greater,
        scaled by the max value in the images.  Otherwise just make sure it
        smooths enough to avoid pixelization effects.

    output
    ------
    width: number
        sqrt(Area)/pi where Area is,

            nim=image.image.max()
            arg =  (nim-thresh)/smooth
            vals = 0.5*( 1 + erf(arg) )
            area = vals.sum()
            width = 2*sqrt(area/pi)
    """
    from numpy import array, sqrt, zeros, pi, where
    from scipy.special import erf


    nim = image.copy()
    maxval=image.max()
    nim *= (1./maxval)

    arg = (nim-thresh)/smooth

    vals = 0.5*( 1+erf(arg) )
    area = vals.sum()
    width = 2*sqrt(area/pi)

    return width

def get_dims_cen(T, nsigma_render):
    """
    Based on T, get the required dimensions and a center
    """
    sigma=numpy.sqrt(T/2.)
    dims = [2.*sigma*nsigma_render]*2
    cen = [(dims[0]-1.)/2.]*2

    return dims, cen


