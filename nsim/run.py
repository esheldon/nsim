from __future__ import print_function

import fitsio
from . import sime
from . import fitters
from . import moments
from . import util

def profile_sim(seed,sim_conf,run_conf,npairs,output_file):
    # don't want to see the JIT
    import cProfile
    import pstats

    import ngmix

    cProfile.runctx('go(seed,sim_conf,run_conf,npairs,output_file)',
                    globals(),locals(),
                    'profile_stats')
    
    p = pstats.Stats('profile_stats')
    p.sort_stats('time').print_stats()

def go(seed,
       sim_conf,
       run_conf,
       npairs,
       output_file,
       make_plots=False,
       write_local=False):


    if seed is None:
        raise ValueError("send a seed")

    sim = sime.Sim(sim_conf, seed)

    plot_base=output_file.replace('.fits','')

    kw=dict(make_plots=make_plots,
            plot_base=plot_base)

    ftype=run_conf['fitter']
    if ftype=='max':
        fitclass=fitters.MaxFitter

    elif ftype=='galsim-max':
        if run_conf['fit_model']=='spergel':
            fitclass=fitters.SpergelFitter
        else:
            fitclass=fitters.GalsimFitter

    elif ftype=='galsim-metacal':
        if run_conf['fit_model']=='spergel':
            fitclass=fitters.SpergelMetacalFitter
        else:
            fitclass=fitters.GalsimMetacalFitter


    elif ftype=='isample-gaussmom':
        fitclass=fitters.ISampleGaussMom

    elif ftype=='psample-gaussmom':
        fitclass=fitters.PSampleGaussMom


    elif ftype=='isample-gaussmom-metacal':
        fitclass=fitters.ISampleGaussMomMetacal


    elif ftype in ['metacal-max','metacal-new']:
        fitclass=fitters.MaxMetacalFitter

    elif ftype == 'metacal-metamom':
        fitclass=fitters.MetacalMetaMomFitter

    elif ftype=='metacal-old':
        fitclass=fitters.MaxMetacalFitterOld

    elif ftype=='metacal-round-analytic-psf':
        fitclass=fitters.MaxMetacalRoundAnalyticPSFFitter


    elif ftype=='nuller-gauss2d':
        fitclass=fitters.NullGauss2Fitter

    elif ftype=='deconv':
        fitclass=fitters.Deconvolver

    elif ftype=='metacal-deconv':
        fitclass=fitters.DeconvMetacalFitter

    elif ftype=='metacal-moments':
        fitclass=fitters.MetacalMomentsOld


    elif ftype=='metacal-momnull-gauss':
        fitclass=moments.MetacalNullGaussFitter


    elif ftype=='metacal-moments-fixed':
        fitclass=moments.MetacalMomentsFixed

    elif ftype=='am':
        fitclass=moments.AMFitter

    elif ftype=='metacal-moments-am':
        # different am fitter
        fitclass=moments.MetacalMomentsAM

    elif ftype=='metacal-moments-am-fixed':
        # different am fitter
        fitclass=moments.MetacalMomentsAMFixed

    elif ftype=='metacal-gaussk':
        # gauss fitter in k space
        fitclass=moments.MetacalGaussK



    elif ftype=='metacal-moments-am-mofsub':
        fitclass=moments.MetacalMomentsAMMOFSub


    elif ftype=='metacal-moments-deweight':
        # different am fitter
        fitclass=moments.MetacalMomentsDeweight



    elif ftype=='metacal-max-simn':
        fitclass=fitters.MaxMetacalSimnFitter
    elif ftype=='metacal-max-subn':
        fitclass=fitters.MaxMetacalSubnFitter

    elif ftype=='metacal-max-fixR':
        fitclass=fitters.MaxMetacalFixRFitter

    elif ftype=='metacal-max-detrend':
        fitclass=fitters.MaxMetacalDetrendFitter


    elif ftype=='metacal-max-degrade':
        if simtype=='galsim':
            fitclass=fitters.MaxMetacalFitterDegradeGS
        else:
            fitclass=fitters.MaxMetacalFitterDegrade

    elif ftype=='ncal':
        fitclass=fitters.NCalFitter

    elif ftype=='pcal':
        fitclass=fitters.PostcalFitter
    elif ftype=='pcal-simn':
        fitclass=fitters.PostcalSimnFitter
    elif ftype=='pcal-simp':
        fitclass=fitters.PostcalSimpFitter

    elif ftype=='pcal-sim-shearp':
        fitclass=fitters.PostcalSimShearpFitter


    elif ftype=='ppmcal':
        fitclass=fitters.PPMetacalFitter

    elif ftype=='metanoise-max':
        if simtype=='galsim':
            fitclass=fitters.MaxMetanoiseFitter
        else:
            raise RuntimeError("don't have metanoise for ngmix sim yet")

    elif ftype=='metacal-metanoise-mom':
        if simtype=='galsim':
            fitclass=fitters.MomMetacalMetanoiseFitter
        else:
            raise RuntimeError("don't have metanoise for ngmix sim yet")



    elif ftype=='metacal-mom':
        if simtype=='galsim':
            fitclass=fitters.MomMetacalFitter
        else:
            raise RuntimeError("don't have mom metacal for ngmix sim yet")

    elif ftype=='metacal-admom':
        fitclass=fitters.AMMetacalFitter


    elif ftype=='metacal-isample':
        fitclass=fitters.ISampleMetacalFitter

    elif ftype=='metacal-isample-nearest':
        fitclass=fitters.ISampleMetacalFitterNearest

    elif ftype=='metacal-psamp':
        fitclass=fitters.PSampleMetacalFitter


    elif ftype=='metacal-em':
        fitclass=fitters.EMMetacalFitter


    elif ftype=='shear-null-postpsf':
        fitclass=fitters.ShearNullFitterPostpsf
    elif ftype=='shear-null-prepsf':
        fitclass=fitters.ShearNullFitterPrepsf

    elif ftype=='kmom-metacal':
        fitclass=fitters.KMomMetacalFitterPost
    elif ftype=='kmom-metacal-pre':
        fitclass=fitters.KMomMetacalFitter
    else:
        raise RuntimeError("bad fitter: '%s'" % ftype)

    fitter=fitclass(sim, run_conf, npairs, **kw)
    fitter.go()

    data=fitter.get_data()

    if write_local:
        success=util.write_fits(output_file, data)
    else:
        print("writing:",output_file)
        fitsio.write(output_file, data, clobber=True)
