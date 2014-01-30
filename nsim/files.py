from sys import stderr
import os
from os.path import join as path_join

import yaml
import numpy


def get_s2n_nrepeat(s2n, fac=0.4):
    """
    Number of repeats.  This is not enough now that I'm using the
    matched s/n

    The 0.4 gives rather noisy results for exp but can run less than a day.
    It gives *very* noisy results for dev, need to increase.
    """
    nrep = round( (fac/( s2n/100. )**2) )
    if nrep < 1:
        nrep = 1
    nrep = int(nrep)
    return nrep


def get_config_dir():
    d=os.environ['NSIM_DIR']
    return path_join(d,'share','nsim_config')

def get_config_file(run):
    d=get_config_dir()
    name='%s.yaml' % run
    return path_join(d, name)

def read_config(run):
    """
    run could be 'name' in sim
    """
    import yaml
    f=get_config_file(run)
    with open(f) as fobj:
        c=yaml.load(fobj)
    if 'run' in c:
        n='run'
    else:
        n='name'
    if c[n] != run:
        raise ValueError("%s in config does not match "
                         "itself: '%s' instead of '%s'" % (n,c[n],run))
    return c

def get_default_fs():
    if os.environ.get('NSIM_FS')=='hdfs':
        fs='hdfs'
    else:
        fs='nfs'
    return fs

def get_simdir(fs=None):
    """
    note old style name for compatibility
    """
    if fs=='hdfs':
        dir=os.environ['SHAPESIM_HDFS_DIR']
    else:
        dir=os.environ['SHAPESIM_DIR']

    return dir

def get_run_dir(run, fs=None):
    dir=get_simdir(fs=fs)
    return path_join(dir,run)


def get_condor_dir(run):
    dir=get_run_dir(run)
    dir=path_join(dir, 'condor')
    return dir

def get_condor_job_url(run, filenum=None, missing=False):
    d=get_condor_dir(run)
    if filenum is not None:
        end = '-%03d' % filenum
    else:
        end=''
    if missing:
        end='%s-missing' % end

    fname='{run}{end}.condor'.format(run=run,end=end)
    return path_join(d,fname)

def get_condor_master_url(run):
    d=get_condor_dir(run)
    return path_join(d,'%s.sh' % run)


def get_output_dir(run, sub=None, fs=None):
    dir=get_run_dir(run, fs=fs)
    dir=path_join(dir, 'outputs')
    if sub:
        dir = path_join(dir, sub)
    return dir

def get_output_url(run, is2, ie, itrial=None, fs=None, ext='fits'):
    """

    is2 and ie are the index in the list of s2 and ellip vals for a given run.
    They should be within [0,nums2) and [0,nume)

    Note ie might actually be is2n

    """
    sub=None
    if itrial is not None:
        sub='bytrial'
    dir=get_output_dir(run, sub=sub, fs=fs)
    f='%s-%03i-%03i' % (run,is2,ie)
    if itrial is not None:
        if itrial == '*':
            f += '-*'
        else:
            f += '-%05i' % itrial
    f += '.%s' % ext
    return path_join(dir, f)


def get_averaged_url(run, is2, fs=None, ext='fits'):
    """
    All the trials are averaged in a given s2 bin, and all
    ellip/s2n bins are in a single struct.
    """

    dir=get_output_dir(run, fs=fs)
    f='%s-%03i-avg' % (run,is2)
    f = '%s.%s' % (f,ext)
    return path_join(dir, f)

s2n_ref_bdfg=[ 10, 15, 23, 35, 53, 81, 123, 187, 285, 433, 658, 1000 ]
npair_ref_bdfg=[1240000, 1240000,  711574,  363878,  164300,
               71517,  31589,  13733,   6045,
               3038,   3038,   3038]
err_ref_bdfg=[ 8.64686523781e-05,6.36021322793e-05,5.78428340247e-05,5.44774877257e-05,5.41502412604e-05,5.39701794399e-05,5.35799931151e-05,5.38147319202e-05,5.35173339764e-05,4.97357492734e-05,3.24657802612e-05,2.14956557248e-05]


# from a BA13 prior run, exp galaxy
# /astro/u/esheldon/lensing/shapesim/cbafit-geg02r07/outputs/cbafit-geg02r07-000-avg.rec
s2n_ref_eg_sr2=numpy.array([ 10, 15, 23, 35, 53, 81, 123, 187, 285, 433, 658, 1000 ],dtype='f8')
err_ref_eg_sr2=numpy.array([5.86094527e-05,   4.24659823e-05,   4.15026596e-05,
                            4.10366038e-05,   4.17969006e-05,   4.19735634e-05,
                            4.19306727e-05,   4.19540700e-05,   4.14286130e-05,
                            3.82807253e-05,   2.54193529e-05,   1.65443563e-05])

npair_ref_eg_sr2=numpy.array([4844000, 4844000, 2351416, 1079866,  459928,  196424,   86000,
                              37368,   16320,    8300,    8300,    8300])


# dg at different sigma ratios
s2n_ref_dg_sr2=numpy.array([ 10, 15, 23, 35, 53, 81, 123, 187, 285, 433, 658, 1000 ],dtype='f8')

err_ref_dg_sr2=numpy.array([5.40373117e-05,   4.38981156e-05,   4.29005342e-05,
                            4.20680925e-05,   4.26229933e-05,   4.28001815e-05,
                            4.26874563e-05,   4.26157464e-05,   4.22199810e-05,
                            3.94608855e-05,   2.59115628e-05,   1.70170656e-05])

npair_ref_dg_sr2=numpy.array([12200000, 12200000,  7000848,  3579968,  1616500,   703696,
                              310856,   135176,    59292,    29768,    29768,    29768])

# gg at different sigma ratios
s2n_ref_gg_sr2=numpy.array([ 10, 15, 23, 35, 53, 81, 123, 187, 285, 433, 658, 1000 ],dtype='f8')
err_ref_gg_sr2=numpy.array([4.76479645748e-05, 3.34475296565e-05, 3.21185913536e-05, 3.15101059311e-05,
                            3.19010304786e-05, 3.19790601662e-05, 3.18597748978e-05, 3.18364047527e-05,
                            3.18242541076e-05, 2.91863314321e-05, 1.91562340848e-05, 1.26330439144e-05])
npair_ref_gg_sr2 = numpy.array([2412000, 2412000, 1169971,  537207,  230174,   98277,   42930,
                                18600,    8140,    4152,    4152,    4152])

s2n_ref_gg_sr1=numpy.array([ 10, 15, 23, 35, 53, 81, 123, 187, 285, 433, 658, 1000 ],dtype='f8')
err_ref_gg_sr1=numpy.array([  3.52809767e-05,   3.53038557e-05,   3.53300171e-05,
                            3.53451981e-05,   3.52917768e-05,   3.52483577e-05,
                            3.52699247e-05,   3.50577662e-05,   3.49212158e-05,
                            3.45652208e-05,   3.32411166e-05,   3.13146698e-05])

npair_ref_gg_sr1 = numpy.array([5424000, 2673128, 1195540,  528388,  232328,   99450,   43200,
                                18648,    8170,    3504,    1512,     656])



def get_npair_by_noise(s2n, desired_err, run, sigratio):
    """
    given the desired final error, determine the required number of pairs
    """

    if '-eg' in run:
        if sigratio >= 1.99:
            npairii = numpy.interp([s2n], s2n_ref_eg_sr2, npair_ref_eg_sr2)
            errii = numpy.interp([s2n], s2n_ref_eg_sr2, err_ref_eg_sr2)
        else:
            raise ValueError("implement other sratio")
    elif 'bdfg' in run:
        # for bdf make sure you set sec_per_pair
        npairii = numpy.interp([s2n], s2n_ref_bdfg, npair_ref_bdfg)
        errii = numpy.interp([s2n], s2n_ref_bdfg, err_ref_bdfg)
    elif '-dg' in run:
        if sigratio >= 1.99:
            npairii = numpy.interp([s2n], s2n_ref_dg_sr2, npair_ref_dg_sr2)
            errii = numpy.interp([s2n], s2n_ref_dg_sr2, err_ref_dg_sr2)
        else:
            raise ValueError("implement other sratio")
    elif '-gg' in run:
        if sigratio >= 1.99:
            npairii = numpy.interp([s2n], s2n_ref_gg_sr2, npair_ref_gg_sr2)
            errii = numpy.interp([s2n], s2n_ref_gg_sr2, err_ref_gg_sr2)
        elif sigratio >= 1.41:
            raise ValueError("implement 1.4")
        else:
            npairii = numpy.interp([s2n], s2n_ref_gg_sr1, npair_ref_gg_sr1)
            errii = numpy.interp([s2n], s2n_ref_gg_sr1, err_ref_gg_sr1)
    else:
        raise ValueError("support runs of type '%s'" % run)

    # desired_err = errii*sqrt(npairii/npair)
    # thus npair = npairii*(errii/desired_err)^2
    npair = npairii*(errii/desired_err)**2

    return npair[0]
    

def get_npair_nsplit_by_noise(c, is2n, npair_min=None):
    """
    Get the npair/nsplit given the input config and s2n index
    """
    from math import ceil

    # this is the requirement from measurement error
    s2n = c['s2n_vals'][is2n]
    sigratio = numpy.sqrt( c['simc']['obj_T_mean']/c['simc']['psf_T'] )
    npair_tot = get_npair_by_noise(s2n, c['desired_err'],c['run'], sigratio)

    npair_shapenoise = 0
    ring=c.get('ring',True)
    if not ring:
        # add in shape noise term for non-ring test
        # err = 0.16/sqrt(ngal)
        ngal = (0.16/c['desired_err'])**2
        npair_shapenoise = int(ngal/2.)


    #print 'meas noise npair:',npair_tot
    #print 'shapenoise npair:',npair_shapenoise 
    npair_tot += npair_shapenoise

    if npair_min is not None and npair_tot < npair_min:
        npair_tot = npair_min

    #print 'desired_err:',c['desired_err']
    #print 'npair_tot:',npair_tot

    # to keep equal time, normalize to zeroth
    nsplit0 = c['nsplit0']

    if is2n==0:
        nsplit=nsplit0
    else:
        npair_tot0 = get_npair_by_noise(c['s2n_vals'][0], c['desired_err'],c['run'],sigratio)
        npair_tot0 += npair_shapenoise
        nsplit = int( ceil( nsplit0*float(npair_tot)/npair_tot0 ))

    npair_per = int(ceil(npair_tot/float(nsplit)))

    return npair_per, nsplit

def get_npair_nsplit(c, is2n, npair_min=None):
    """
    Get number of pairs per split and number of splits

    For equal_time, we take number per split from is2n==0
    """
    if 'desired_err' in c:
        return get_npair_nsplit_by_noise(c, is2n, npair_min=npair_min)
    else:
        s2n = c['s2n_vals'][is2n]

        npair = get_s2n_nrepeat(s2n, fac=c['s2n_fac'])
        if npair < c['min_npair']:
            npair = c['min_npair']
        nsplit=c['nsplit']

        return npair,nsplit

