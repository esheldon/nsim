from __future__ import print_function
import galsim
import ngmix
from . import pdfs

def get_psf_maker(config, rng):
    if config['model'] == 'multi-component':
        maker=MultiComponentPSF(config, rng)
    else:
        maker=PSFSampler(config, rng)

    return maker

class PSFSampler(dict):
    """
    This basic sampler does gauss and Moffat

    psf subset of the config.

    # Gaussian with fixed properties
    model: 'gauss'
    r50:
        1.5
    shape:
        [0.0, 0.025]


    # Moffat with variable properties
    model: 'moffat'
    beta: 3.5
    r50:
        type: "lognormal"
        mean: 2.0
        sigma: 0.2
    shape:
        type: "normal2d"
        cen: [0.0, 0.025]
        sigma: [0.01, 0.01]
    """
    def __init__(self, config, rng):
        self.update(config)
        self.rng=rng
        self._set_pdfs()

    def __call__(self):

        model=self['model']

        kw = self._get_keywords()

        if model=='moffat':
            kw['beta'] = self['beta']
            cls = galsim.Moffat
        elif model=='gauss':
            cls = galsim.Gaussian
        else:
            raise ValueError("bad psf model: '%s'" % model)

        psf = cls(**kw)

        g1, g2 = self._get_shape()
        psf = psf.shear(g1=g1, g2=g2)

        return psf, kw

    def _get_shape(self):
        if self.shape_pdf is not None:
            try:
                g1, g2 = self.shape_pdf.sample()
            except:
                g1, g2 = self.shape_pdf.sample2d()

            print("    psf (pdf) shape: %g %g" % (g1, g2))

        else:
            g1=self['shape'][0]
            g2=self['shape'][1]

        return g1, g2

    def _get_keywords(self):

        kw={}

        if 'r50' in self:
            if self.r50_pdf is not None:
                r50 = self.r50_pdf.sample()
                print("    psf r50: %g" % r50)

            elif 'r50' in self:
                r50 = self['r50']

            else:
                raise ValueError("r50 value or distribution must be set")

            kw['half_light_radius'] = r50

        elif 'fwhm' in self:
            if self.fwhm_pdf is not None:
                fwhm = self.fwhm_pdf.sample()
                print("    psf fwhm: %g" % fwhm)

            elif 'fwhm' in self:
                fwhm = self['fwhm']

            else:
                raise ValueError("fwhm value or distribution must be set")

            kw['fwhm'] = fwhm

        return kw


    def _set_pdfs(self):

        self._set_shape_pdf()
        self._set_size_pdf()

    def _set_shape_pdf(self):
        self.shape_pdf=None

        shapeconf = self['shape']

        if isinstance(shapeconf,dict):
            assert shapeconf['type']=="normal2d"

            if shapeconf['type']=='normal2d':
                self.shape_pdf=ngmix.priors.SimpleGauss2D(
                    shapeconf['cen'][0],
                    shapeconf['cen'][1],
                    shapeconf['sigma'][0],
                    shapeconf['sigma'][1],
                    rng=self.rng,
                )
        else:
            if len(shapeconf) != 2:
                raise ValueError("for constant psf "
                                 "shapes, length must be 2")

    def _set_size_pdf(self):
        if 'r50' in self:
            self._set_r50_pdf()
        elif 'fwhm' in self:
            self._set_fwhm_pdf()
        else:
            raise ValueError("need r50 or fwhm in psf config")

    def _set_r50_pdf(self):
        self.r50_pdf=None
        r50conf = self['r50']

        if isinstance(r50conf, dict):

            if r50conf['type'] == 'lognormal':
                self.r50_pdf = ngmix.priors.LogNormal(
                    r50conf['mean'],
                    r50conf['sigma'],
                    rng=self.rng,
                )
            else:
                raise ValueError("bad psf r50 pdf "
                                 "type: '%s'" % r50conf['type'])

    def _set_fwhm_pdf(self):
        self.fwhm_pdf=None
        fwhmconf = self['fwhm']

        if isinstance(fwhmconf, dict):

            if fwhmconf['type'] == 'lognormal':
                self.fwhm_pdf = ngmix.priors.LogNormal(
                    fwhmconf['mean'],
                    fwhmconf['sigma'],
                    rng=self.rng,
                )
            else:
                raise ValueError("bad psf fwhm pdf "
                                 "type: '%s'" % fwhmconf['type'])


class MultiComponentPSF(object):
    def __init__(self, config, rng):
        self.rng=rng
        self._set_components(config)

    def __call__(self):

        parts=[]
        for c in self.components:
            tobj, junk = c()
            parts.append(tobj)

        obj = galsim.Add(parts)

        r50=obj.calculateHLR()
        meta = {'r50':r50}
        return obj, meta 

    def _set_components(self, config):
        self.config=config

        pieces=[]
        for component in config['components']:

            pieces.append( PSFSampler(component, self.rng) )

        self.components=pieces

