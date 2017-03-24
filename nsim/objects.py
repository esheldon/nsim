from __future__ import print_function
import galsim
import ngmix
from . import pdfs

def get_object_maker(config, rng):
    model=config['model']
    if model in ['gauss','exp','dev']:
        maker=SimpleObjectMaker(config, rng)
    else:
        raise ValueError("bad model: '%s'" % model)

class SimpleObjectMaker(dict):
    """
    simple models such as exp.  The settable parameters are only

        - size
        - flux
        - g1,g2

    Things like shift within the stamp are done in the image maker,
    shear is also done elsewhere 
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
        return gal, {'r50':r50,'flux':flux,'g1':g1,'g2':g2}


    def _set_pdf(self):
        """
        set distributions for each parameter

        We use joint distributions for flux and size, though these could
        be separable in practice
        """

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
            raise ValueError("only separable pdfs currently supported")


    def _get_g_pdf(self):
        if 'g' in self['obj_model']:
            g_spec=self['obj_model']['g']
            self.g_pdf=ngmix.priors.GPriorBA(
                g_spec['sigma'],
                rng=self.rng,
            )
        else:
            self.g_pdf=None

    def _get_r50_pdf(self):
        r50spec = self['obj_model']['r50']
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

        fluxspec = self['obj_model']['flux']
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
