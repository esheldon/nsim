from __future__ import print_function
import logging
import galsim
import ngmix
from . import pdfs

logger = logging.getLogger(__name__)

def get_object_maker(config, rng, galsim_rng):
    model=config['model']
    if model in ['gauss','exp','dev']:
        maker=SimpleMaker(config, rng)
    elif model in ['bdk','bd']:
        maker=BDKMaker(config, rng, galsim_rng)
    elif model in ['knots']:
        maker=KnotsMaker(config, rng, galsim_rng)
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
        return self._make_object(**kw)

    def _make_object(self, **kw):
        no_central=kw.get('no_central',False)
        no_central=kw.get('no_central',None)
        if no_central is None:
            no_central=self.get('no_central',False)

        if 'nbrs' in self:
            nconf=self['nbrs']

            objs=[]
            if not no_central:
                logger.debug('    adding central')
                central, meta = self._get_one_model(**kw)
                objs.append(central)
            else:
                logger.debug('    not adding central')

            for i in range(nconf['num']):
                nbr, nmeta = self._get_one_model(**kw)
                if 'flux_frac' in nconf:
                    nbr = nbr.withScaledFlux(nconf['flux_frac'])
                else:
                    nbr = nbr.withFlux(nconf['flux'])

                dx, dy = self._get_nbr_shift()
                nbr=nbr.shift(dx=dx, dy=dy)

                objs.append(nbr)

            obj = galsim.Sum(objs)

            if no_central:
                meta=nmeta

        else:
            obj,meta = self._get_one_model(**kw)

        return obj, meta

    def _get_nbr_shift(self):
        if hasattr(self._nbr_shift_pdf,'sample2d'):
            dx, dy = self._nbr_shift_pdf.sample2d()
        else:
            dx, dy = self._nbr_shift_pdf.sample(2)

        return dx, dy

    def _get_one_model(self, **kw):
        g1,g2,r50,flux = self.pdf.sample()

        if 'flux' in kw:
            flux=kw['flux']
        if 'r50' in kw:
            r50=kw['r50']

        fracknots, nknots = self._get_knot_info()

        # if knots were not requested, fracknots defaults to 0.0
        primary_flux = (1.0 - fracknots)*flux

        model=self['model']
        if model=='gauss':
            gal = galsim.Gaussian(flux=primary_flux, half_light_radius=r50)
        elif model=='exp':
            gal = galsim.Exponential(flux=primary_flux, half_light_radius=r50)
        elif model=='dev':
            gal = galsim.DeVaucouleurs(flux=primary_flux, half_light_radius=r50)
        else:
            raise ValueError("bad galaxy model: '%s'" % pars['model'])

        if nknots > 0:
            knot_flux = fracknots*flux

            knots = galsim.RandomWalk(
                npoints=nknots,
                half_light_radius=r50,
                flux=knot_flux,
                rng=self.galsim_rng,
            )

            gal = galsim.Add([gal, knots])

        gal = gal.shear(g1=g1, g2=g2)
        meta={
            'r50':r50,
            'flux':flux,
            'fracknots':fracknots,
            'nknots':nknots,
        }

        return gal, meta

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
        elif 'r50_flux' in self:
            # joint in r50 and flux

            g_pdf = self._get_g_pdf()
            r50_flux_pdf = self._get_joint_r50_flux_pdf()

            self.pdf = pdfs.ShapeJointR50FluxPDF(
                r50_flux_pdf,
                g_pdf=g_pdf,
            )

        else:
            raise ValueError("only separable flux/size/shape "
                             "pdfs currently supported")

        self.fracknots_pdf = self._get_fracknots_pdf()

        if 'nbrs' in self:
            self._nbr_shift_pdf=self._get_nbr_shift_pdf()

    def _get_nbr_shift_pdf(self):
        sconf=self['nbrs']['shift']

        if sconf['type']=='uniform':
            return ngmix.priors.FlatPrior(
                -sconf['radius'], sconf['radius'],
                rng=self.rng,
            )
        elif sconf['type']=='disk':
            return ngmix.priors.ZDisk2D(
                sconf['radius'],
                rng=self.rng,
            )
        else:
            raise ValueError('bad nbr shift pdf: %s' % sconf['type'])

    def _get_joint_r50_flux_pdf(self):
        """
        joint size-flux from the cosmos catalog
        """

        spec=self['r50_flux']
        if spec['type']=='cosmos':
            flux_range = spec.get('flux_range', [2.5, 100.0])
            r50_range = spec.get('r50_range', [0.15, 1.0])


            pdf = pdfs.CosmosR50Flux(
                r50_range,
                flux_range,
            )

        else:
            raise ValueError("bad r50_flux joint "
                             "pdf: '%s'" % spec['type'])

        return pdf


    def _get_g_pdf(self):
        if 'g' not in self:
            return None

        g_spec=self['g']

        if g_spec['type']=="ba":
            g_pdf=ngmix.priors.GPriorBA(
                g_spec['sigma'],
                rng=self.rng,
            )
        elif g_spec['type'] == "gauss":
            g_pdf=ngmix.priors.GPriorGauss(
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
                logger.debug("Reading r50 values from file: %s" % fname)
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

                flux_pdf=pdfs.PowerLaw(index, xmin, xmax)

            else:
                raise ValueError("bad flux pdf type: '%s'" % fluxspec['type'])

        return flux_pdf

    def _get_knot_info(self):
        if self.fracknots_pdf is None:
            fracknots=0.0
            num=0
        else:
            fracknots = self.fracknots_pdf.sample()
            num = self['knots']['num']

        return fracknots, num

    def _get_fracknots_pdf(self):
        if 'knots' not in self:
            kconf=None
        else:
            kconf = self['knots']

        if kconf is None:
            pdf = None
        else:
            bdr = self['knots']['flux_frac']['range']
            pdf = ngmix.priors.FlatPrior(
                bdr[0],
                bdr[1],
                rng=self.rng,
            )

        return pdf

class KnotsMaker(SimpleMaker):
    """
    pure knots
    """
    def __init__(self, config, rng, galsim_rng):
        super(KnotsMaker,self).__init__(config, rng)

        self.galsim_rng=galsim_rng


    def _get_one_model(self, **kw):
        g1,g2,r50,flux = self.pdf.sample()

        if 'flux' in kw:
            flux=kw['flux']
        if 'r50' in kw:
            r50=kw['r50']

        gal = galsim.RandomWalk(
            npoints=self['num'],
            half_light_radius=r50,
            flux=flux,
            rng=self.galsim_rng,
        )

        gal = gal.shear(g1=g1, g2=g2)
        meta={
            'r50':r50,
            'flux':flux,
            'nknots':self['num'],
        }

        return gal, meta

class BDKMaker(SimpleMaker):
    def __init__(self, config, rng, galsim_rng):
        super(BDKMaker,self).__init__(config, rng)

        self.galsim_rng=galsim_rng

    def _set_pdf(self):
        super(BDKMaker,self)._set_pdf()
        self.fracdev_pdf = self._get_fracdev_pdf()

        # component shifts in sky units
        self.shift_pdf = self._get_shift_pdf()

        self.bulge_size_fac_pdf=None
        spec = self.get('bulge_size_fac',None)
        if spec is not None:
            assert spec['type']=='uniform'
            rng=spec['range']
            self.bulge_size_fac_pdf= ngmix.priors.FlatPrior(
                rng[0], rng[1],
                rng=self.rng,
            )

    def _get_component_sizes(self, r50):

        disk_r50=r50

        if self.bulge_size_fac_pdf is not None:
            fac = self.bulge_size_fac_pdf.sample()
        else:
            fac=1

        bulge_r50 = disk_r50*fac
        return disk_r50, bulge_r50

    def _get_one_model(self, **kw):
        g1disk,g2disk,r50,flux = self.pdf.sample()
        g1bulge,g2bulge,junk,junk = self.pdf.sample()

        if 'flux' in kw:
            flux=kw['flux']
        if 'r50' in kw:
            r50=kw['r50']

        disk_r50, bulge_r50 = self._get_component_sizes(r50)

        fracdev = self.fracdev_pdf.sample()

        fracknots, nknots = self._get_knot_info()
        shifts = self._get_component_shifts()

        bulge_flux = flux*fracdev

        disk_flux_total = flux*(1.0 - fracdev)
        disk_flux = (1.0 - fracknots)*disk_flux_total
        knot_flux = fracknots*disk_flux_total

        disk_raw = galsim.Exponential(
            flux=disk_flux,
            half_light_radius=disk_r50,
        )

        if nknots > 0:
            knots = galsim.RandomWalk(
                npoints=nknots,
                half_light_radius=disk_r50,
                flux=knot_flux,
                rng=self.galsim_rng,
            )

            disk = galsim.Add([disk_raw, knots])
        else:
            disk = disk_raw

        # the bulge is always smooth
        bulge = galsim.DeVaucouleurs(
            flux=bulge_flux,
            half_light_radius=bulge_r50,
        )

        # disk and bulge get independent shapes
        disk  = disk.shear(g1=g1disk, g2=g2disk)
        bulge = bulge.shear(g1=g1bulge, g2=g2bulge)

        if shifts is not None:
            #print("shifting disk: ",shifts['disk'])
            #print("shifting bulge:",shifts['bulge'])
            disk  = disk.shift(
                dx=shifts['disk']['dx'],
                dy=shifts['disk']['dy'],
            )
            bulge  = bulge.shift(
                dx=shifts['bulge']['dx'],
                dy=shifts['bulge']['dy'],
            )


        # combine them and shear that
        gal = galsim.Add([disk, bulge])

        meta={
            'r50':r50,
            'flux':flux,
            'fracdev':fracdev,
            'fracknots':fracknots,
            'nknots':nknots,
        }

        return gal, meta

    def _get_fracdev_pdf(self):
        bdr = self['fracdev']['range']
        return ngmix.priors.FlatPrior(bdr[0], bdr[1],
                                      rng=self.rng)

    def _get_component_shifts(self):
        if 'component_offsets' in self:
            dx,dy = self.shift_pdf.sample2d()
            shift_disk = {'dx':dx, 'dy':dy}
            dx,dy = self.shift_pdf.sample2d()
            shift_bulge = {'dx':dx, 'dy':dy}
        else:
            shift_disk = {'dx':0.0, 'dy':0.0}
            shift_bulge = {'dx':0.0, 'dy':0.0}

        return {
            'disk':shift_disk,
            'bulge':shift_bulge,
        }
    def _get_shift_pdf(self):
        if 'component_offsets' in self:
            return ngmix.priors.ZDisk2D(
                self['component_offsets']['radius'],
                rng=self.rng,
            )
        else:
            return None
