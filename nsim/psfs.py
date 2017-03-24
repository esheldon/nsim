from __future__ import print_function
import galsim
import ngmix
from . import pdfs

def get_psf_maker(config, rng):
    if config]['model'] == 'multi-component':
        maker=psfs.MultiComponentPSF(config, rng)
    else:
        maker=psfs.PSFSampler(config, rng)

    return maker

class PSFSampler(object):
    """
    Basic sampler
    """
    def __init__(self, config, rng):
        self.rng=rng
        self.set_pdfs(config)

    def __call__(self):
        config = self.config

        model=config['model']

        flux = config.get('flux',1.0)

        r50, fwhm = self.get_size()

        if model=='moffat':
            psf = galsim.Moffat(
                beta=config['beta'],
                half_light_radius=r50,
                fwhm=fwhm,
                flux=flux,
            )
        elif model=='gauss':
            psf = galsim.Gaussian(
                half_light_radius=r50,
                fwhm=fwhm,
                flux=flux,
            )
        else:
            raise ValueError("bad psf model: '%s'" % model)


        psf_g1, psf_g2 = self.get_shape()
        psf = psf.shear(g1=psf_g1, g2=psf_g2)

        return psf, {'fwhm':fwhm, 'r50':r50}

    def get_shape(self):
        config = self.config

        if self.psf_ellip_pdf is not None:
            try:
                psf_g1, psf_g2 = self.psf_ellip_pdf.sample()
            except:
                psf_g1, psf_g2 = self.psf_ellip_pdf.sample2d()

            print("    psf (pdf) shape: %g %g" % (psf_g1, psf_g2))

        elif 'randomized_orientation' in config:
            ro=config['randomized_orientation']
            if ro["dist"]=="uniform":
                #angle = numpy.random.random()*2*numpy.pi
                angle = self.rng.uniform()*2*numpy.pi
                psf_shape = ngmix.Shape(ro['magnitude'], 0.0)
                psf_shape.rotate(angle)
                psf_g1 = psf_shape.g1
                psf_g2 = psf_shape.g2
                print("    psf rand orient. shape: %g %g" % (psf_g1, psf_g2))
            else:
                raise ValueError("only uniform randomized psf orientation for now")
        else:
            psf_g1=config['shape'][0]
            psf_g2=config['shape'][1]

        return psf_g1, psf_g2

    def get_size(self):
        config=self.config

        r50=None
        fwhm=None

        if self.psf_r50_pdf is not None:
            r50 = self.psf_r50_pdf.sample()
            print("    psf r50: %g" % r50)

        elif self.psf_fwhm_pdf is not None:
            fwhm = self.psf_fwhm_pdf.sample()
            print("    psf fwhm: %g" % fwhm)

        elif 'r50' in config:
            r50 = config['r50']

        elif 'fwhm' in config:
            fwhm = config['fwhm']

        return r50, fwhm


    def set_pdfs(self, config):

        self.config=config

        self.psf_r50_pdf=None
        self.psf_fwhm_pdf=None

        if 'fwhm' in config:
            if isinstance(config['fwhm'], dict):
                type=config['fwhm']['type']
                if type=='discrete-pdf':

                    fname=os.path.expandvars( config['fwhm']['file'] )
                    print("Reading fwhm values from file:",fname)
                    vals=numpy.fromfile(fname, sep='\n')
                    psf_fwhm_pdf=psfs.DiscreteSampler(
                        vals,
                        rng=self.rng,
                    )

                elif type=='lognormal':
                    psf_fwhm_pdf=ngmix.priors.LogNormal(
                        config['fwhm']['mean'],
                        config['fwhm']['sigma'],
                        rng=self.rng,
                    )
                else:
                    raise ValueError("bad fwhm type: '%s'" % type)

                if 'range' in config['fwhm']:
                    bounds=config['fwhm']['range']
                    print("   bounding psf fwhm to:",bounds)
                    psf_fwhm_pdf = ngmix.priors.Bounded1D(psf_fwhm_pdf,bounds)

                self.psf_fwhm_pdf  = psf_fwhm_pdf 
        else:
            if isinstance(config['r50'], dict):
                r50pdf = config['r50']
                assert r50pdf['type']=="lognormal","r50 pdf log normal for now"

                self.psf_r50_pdf = ngmix.priors.LogNormal(
                    r50pdf['mean'],
                    r50pdf['sigma'],
                    rng=self.rng,
                )

        if 'shapenoise' in config:
            self.psf_ellip_pdf = ngmix.priors.GPriorBA(config['shapenoise'])

        elif isinstance(config['shape'],dict):
            ppdf=config['shape']
            assert ppdf['type']=="normal2d"
            self.psf_ellip_pdf=ngmix.priors.SimpleGauss2D(
                ppdf['cen'][0],
                ppdf['cen'][1],
                ppdf['sigma'][0],
                ppdf['sigma'][1],
                rng=self.rng,
            )
        else:
            self.psf_ellip_pdf=None


class MultiComponentPSF(object):
    def __init__(self, config, rng):
        self.rng=rng
        self.set_components(config)

    def __call__(self):

        parts=[]
        for c in self.components:
            tobj, junk = c()
            parts.append(tobj)

        obj = galsim.Add(parts)

        r50=obj.calculateHLR()
        fwhm=obj.calculateFWHM()
        meta = {'r50':r50, 'fwhm':fwhm}
        return obj, meta 

    def set_components(self, config):
        self.config=config

        pieces=[]
        for component in config['components']:

            pieces.append( PSFSampler(component, self.rng) )

        self.components=pieces

