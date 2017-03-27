from __future__ import print_function
import galsim
import ngmix
from . import pdfs
from .shearpdf import get_shear_pdf

def get_object_maker(config, rng):
    model=config['model']
    if model in ['gauss','exp','dev']:
        maker=SimpleObjectMaker(config, rng)
    else:
        raise ValueError("bad model: '%s'" % model)

class SimpleObjectMaker(dict):
    """
    object section of config

    Things like shift within the stamp are done in the image maker,
    shear is also done elsewhere 

    model: 'exp'
    flux:
        type: "lognormal"
        mean: 100.0
        sigma: 30.0

    r50:
        type: "lognormal"
        mean: 2.0
        sigma: 1.0

    g:
        type: "ba"
        sigma: 0.2

    shear:
        shears: [0.02, 0.00]

    """
    def __init__(self, config, rng):
        self.update(config)
        self.rng=rng

        self._set_pdf()

    def __call__(self):
        g1,g2,r50,flux = self.pdf.sample()

        if pars['model']=='gauss':
            gal = galsim.Gaussian(flux=flux, half_light_radius=r50)
        elif pars['model']=='exp':
            gal = galsim.Exponential(flux=flux, half_light_radius=r50)
        elif pars['model']=='dev':
            gal = galsim.DeVaucouleurs(flux=flux, half_light_radius=r50)
        else:
            raise ValueError("bad galaxy model: '%s'" % pars['model'])

        gal = gal.shear(g1=g1, g2=g2)

        if self.shear_pdf is not None:
            shear, shindex = self.shear_pdf.get_shear(self.rng)
            gal = gal.shear(g1=shear.g1, g2=shear.g2)
        else:
            shindex=0

        return gal, {'r50':r50,'flux':flux,'g1':g1,'g2':g2,'shear_index':shindex}


    def _set_pdf(self):
        """
        set distributions for each parameter

        We use joint distributions for flux and size, though these could
        be separable in practice
        """

        
        if 'shear' in self:
            self.shear_pdf = get_shear_pdf(self['shear'])
        else:
            self.shear_pdf = None

        if 'flux' in self and 'r50' in self:
            # the pdfs are separate

            g_pdf = self._get_g_pdf()
            r50_pdf = self._get_r50_pdf()
            flux_pdf = self._get_flux_pdf()

            self.pdf = pdfs.SeparableShapeR50FluxPDF(
                r50_pdf,
                flux_pdf,
                g_pdf=g_pdf,
            )
        else:
            raise ValueError("only separable flux/size/shape pdfs currently supported")


    def _get_g_pdf(self):
        g_spec=self.get('g',None)
        g_pdf=None

        if g_spec is not None:

            if g_spec['type']=="ba":
                g_pdf=ngmix.priors.GPriorBA(
                    g_spec['sigma'],
                    rng=self.rng,
                )
            else:
                raise ValueError("bad g type: '%s'" % g_spec['type'])

        return g_pdf

    def _get_r50_pdf(self):

        r50spec = self['r50']

        if not isinstance(r50spec,dict):
            r50_pdf=DiscreteSampler([r50spec], rng=self.rng)
        else:

            if r50spec['type']=='uniform':
                r50_r = r50spec['range']
                r50_pdf=ngmix.priors.FlatPrior(
                    r50_r[0],
                    r50_r[1],
                    rng=self.rng,
                )
            elif r50spec['type']=='lognormal':
                r50_pdf=ngmix.priors.LogNormal(
                    r50spec['mean'],
                    r50spec['sigma'],
                    rng=self.rng,
                )
            elif r50spec['type']=='discrete-pdf':
                fname=os.path.expandvars( r50spec['file'] )
                print("Reading r50 values from file:",fname)
                vals=fitsio.read(fname)
                r50_pdf=DiscreteSampler(vals, rng=self.rng)

            else:
                raise ValueError("bad r50 pdf type: '%s'" % r50spec['type'])

        return r50_pdf

    def _get_flux_pdf(self):

        fluxspec = self['flux']
        if not isinstance(fluxspec,dict):
            flux_pdf=DiscreteSampler([fluxspec], rng=self.rng)
        else:

            if fluxspec['type']=='uniform':
                flux_r = fluxspec['range']
                flux_pdf=ngmix.priors.FlatPrior(
                    flux_r[0], flux_r[1],
                    rng=self.rng,
                )
            elif fluxspec['type']=='lognormal':
                flux_pdf=ngmix.priors.LogNormal(
                    fluxspec['mean'],
                    fluxspec['sigma'],
                    rng=self.rng,
                )
            elif fluxspec['type']=='gmixnd':
                flux_pdf=load_gmixnd(fluxspec,rng=self.rng)

            elif fluxspec['type']=='powerlaw':

                index=fluxspec['index']
                xmin=fluxspec['min']
                xmax=fluxspec['max']

                flux_pdf=PowerLaw(index, xmin, xmax)

            else:
                raise ValueError("bad flux pdf type: '%s'" % fluxspec['type'])

        return flux_pdf
