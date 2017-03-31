from __future__ import print_function
import galsim
import ngmix
from . import pdfs

def get_object_maker(config, rng, galsim_rng):
    model=config['model']
    if model in ['gauss','exp','dev']:
        maker=SimpleMaker(config, rng)
    elif model =='bdk':
        maker=BDKMaker(config, rng, galsim_rng)
    else:
        raise ValueError("bad model: '%s'" % model)

    return maker

class SimpleMaker(dict):
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
    """
    def __init__(self, config, rng):
        self.update(config)
        self.rng=rng

        self._set_pdf()

    def __call__(self, **kw):

        gal, meta = self._make_object(**kw)

        return gal, meta

    def _make_object(self, **kw):
        g1,g2,r50,flux = self.pdf.sample()

        if 'flux' in kw:
            flux=kw['flux']
        if 'r50' in kw:
            r50=kw['r50']

        model=self['model']
        if model=='gauss':
            gal = galsim.Gaussian(flux=flux, half_light_radius=r50)
        elif model=='exp':
            gal = galsim.Exponential(flux=flux, half_light_radius=r50)
        elif model=='dev':
            gal = galsim.DeVaucouleurs(flux=flux, half_light_radius=r50)
        else:
            raise ValueError("bad galaxy model: '%s'" % pars['model'])

        gal = gal.shear(g1=g1, g2=g2)
        meta={
            'r50':r50,
            'flux':flux,
        }

        return gal, meta



    def _set_pdf(self):
        """
        set distributions for each parameter

        We use joint distributions for flux and size, though these could
        be separable in practice
        """

        

        if 'flux' in self and 'r50' in self and 'g' in self:
            # the pdfs are separate

            g_pdf = self._get_g_pdf()
            r50_pdf = self._get_r50_pdf()
            flux_pdf = self._get_flux_pdf()

            self.pdf = pdfs.SeparableShapeR50FluxPDF(
                r50_pdf,
                flux_pdf,
                g_pdf=g_pdf,
            )
        elif 'r50_flux' in self:
            # joint in r50 and flux

            g_pdf = self._get_g_pdf()
            r50_flux_pdf = self._get_joint_r50_flux_pdf()

            self.pdf = pdfs.ShapeJointR50FluxPDF(
                r50_pdf,
                flux_pdf,
                g_pdf=g_pdf,
            )

        else:
            raise ValueError("only separable flux/size/shape "
                             "pdfs currently supported")

    def _get_joint_r50_flux_pdf(self):
        """
        joint size-flux from the cosmos catalog
        """

        spec=self['r50_flux']
        if spec['type']=='cosmos':
            pdf = pdfs.CosmosR50Flux(
                spec['r50_range'],
                spec['flux_range'],
            )
        else:
            raise ValueError("bad r50_flux joint "
                             "pdf: '%s'" % spec['type'])

        return pdf


    def _get_g_pdf(self):
        if 'g' not in self:
            raise ValueError("no g spec found in config")

        g_spec=self['g']

        if g_spec['type']=="ba":
            g_pdf=ngmix.priors.GPriorBA(
                g_spec['sigma'],
                rng=self.rng,
            )
        else:
            raise ValueError("bad g type: '%s'" % g_spec['type'])

        return g_pdf

    def _get_r50_pdf(self):

        if 'r50' not in self:
            raise ValueError("no r50 spec found in config")

        r50spec = self['r50']

        if not isinstance(r50spec,dict):
            r50_pdf=pdfs.DiscreteSampler([r50spec], rng=self.rng)
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
                r50_pdf=pdfs.DiscreteSampler(vals, rng=self.rng)

            else:
                raise ValueError("bad r50 pdf type: '%s'" % r50spec['type'])

        return r50_pdf

    def _get_flux_pdf(self):

        if 'flux' not in self:
            raise ValueError("no flux spec found in config")

        fluxspec = self['flux']
        if not isinstance(fluxspec,dict):
            flux_pdf=pdfs.DiscreteSampler([fluxspec], rng=self.rng)
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

class BDKMaker(SimpleMaker):
    def __init__(self, config, rng, galsim_rng):
        super(BDKMaker,self).__init__(config, rng)

        self.galsim_rng=galsim_rng

    def _set_pdf(self):
        super(BDKMaker,self)._set_pdf()

        self.fracdev_pdf = self._get_fracdev_pdf()
        self.fracknots_pdf = self._get_fracknots_pdf()


    def _make_object(self, **kw):
        g1disk,g2disk,r50,flux = self.pdf.sample()
        g1bulge,g2bulge,junk,junk = self.pdf.sample()

        if 'flux' in kw:
            flux=kw['flux']
        if 'r50' in kw:
            r50=kw['r50']


        fracdev = self.fracdev_pdf.sample()
        fracknots = self.fracknots_pdf.sample()

        bulge_flux = flux*fracdev

        disk_flux_total = flux*(1.0 - fracdev)
        disk_flux = (1.0 - fracknots)*disk_flux_total
        knot_flux = fracknots*disk_flux_total

        disk_raw = galsim.Exponential(
            flux=disk_flux,
            half_light_radius=r50,
        )

        knots = galsim.RandomWalk(
            npoints=self['knots']['num'],
            half_light_radius=r50,
            flux=knot_flux,
            rng=self.galsim_rng,
        )

        disk = galsim.Add([disk_raw, knots])

        # the bulge is always smooth
        bulge = galsim.DeVaucouleurs(flux=bulge_flux, half_light_radius=r50)

        # disk and bulge get independent shapes
        disk  = disk.shear(g1=g1disk, g2=g2disk)
        bulge = bulge.shear(g1=g1bulge, g2=g2bulge)

        # combine them and shear that
        gal = galsim.Add([disk, bulge])

        meta={
            'r50':r50,
            'flux':flux,
            'fracdev':fracdev,
            'fracknots':fracknots,
        }

        return gal, meta

    def _get_fracdev_pdf(self):
        bdr = self['fracdev']['range']
        return ngmix.priors.FlatPrior(bdr[0], bdr[1],
                                      rng=self.rng)

    def _get_fracknots_pdf(self):
        bdr = self['knots']['flux_frac']['range']
        return ngmix.priors.FlatPrior(bdr[0], bdr[1],
                                      rng=self.rng)


