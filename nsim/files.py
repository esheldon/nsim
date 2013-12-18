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
    if docum:
        f += '-cum'
    f = '%s.s' % (f,ext)
    return path_join(dir, f)



# from a BA13 prior run, exp galaxy
# /astro/u/esheldon/lensing/shapesim/cbafit-geg02r07/outputs/cbafit-geg02r07-000-avg.rec
s2n_ref_geg=[15,  20,  25,  30,  40,  50,  60,  70,  80,  90,  100,  120,  140,  160,  180,  200,  250,  300,  350,  400]

err_ref_geg=[6.04966948076, 6.27720009086, 6.39230479966, 6.46027631805, 6.31490331065, 5.07534788684, 4.24058192767, 3.63903754734, 3.18832441297, 2.83577943389, 2.55346528496, 2.12927717929, 1.82605936882, 1.5983608504, 1.42135553436, 1.27964356694, 1.02561949518, 0.857061929443, 0.737296676262, 0.648170473541]
npair_ref_geg=[1333000,  750000,  480000,  333000,  200000,  200000,  200000,  200000,  200000,  200000,  200000,  200000,  200000,  200000,  200000,  200000,  200000,  200000,  200000,  200000]

s2n_ref_geg=numpy.array(s2n_ref_geg,dtype='f8')
err_ref_geg=numpy.array(err_ref_geg)*1.e-5
npair_ref_geg=numpy.array(npair_ref_geg)


s2n_ref_deg=numpy.array([15,  20,  25,  30,  40,  50,  60,  70,  80,  90,  100,  120,  140,  160,  180,  200,  250,  300,  350,  400],dtype='f8')


err_ref_deg=numpy.array([  5.72424798e-05,   6.10264401e-05,   6.32783893e-05,
                         6.47026643e-05,   6.63254272e-05,   6.71194744e-05,
                         6.75744094e-05,   6.78835605e-05,   6.80540579e-05,
                         6.81622436e-05,   6.82656513e-05,   6.84032611e-05,
                         6.84544518e-05,   6.84826729e-05,   6.85100035e-05,
                         6.85641846e-05,   6.84825446e-05,   6.83578669e-05,
                         6.82696058e-05,   6.80481557e-05])
npair_ref_deg=numpy.array([3768000, 2280984, 1514069, 1072740,  615414,  397620,  277595,
                           204484,  156980,  124212,  100687,   69975,   51480,   39438,
                           31178,   25272,   16236,   11336,    8397,    6489])

def get_npair_by_noise(s2n, desired_err, run):
    """
    given the desired final error, determine the required number of pairs
    """

    if 'geg' in run or '-eg' in run:
        npairii = numpy.interp([s2n], s2n_ref_geg, npair_ref_geg)
        errii = numpy.interp([s2n], s2n_ref_geg, err_ref_geg)
    elif 'deg' in run or '-dg' in run:
        npairii = numpy.interp([s2n], s2n_ref_deg, npair_ref_deg)
        errii = numpy.interp([s2n], s2n_ref_deg, err_ref_deg)

    # desired_err = errii*sqrt(npairii/npair)
    # thus npair = npairii*(errii/desired_err)^2
    npair = npairii*(errii/desired_err)**2

    return npair[0]
    

def get_npair_nsplit_by_noise(c, is2n):
    """
    Get the npair/nsplit given the input config and s2n index
    """
    from math import ceil
    s2n = c['s2n_vals'][is2n]
    npair_tot = get_npair_by_noise(s2n, c['desired_err'],c['run'])
    #print 'desired_err:',c['desired_err']
    #print 'npair_tot:',npair_tot

    # to keep equal time, normalize to zeroth
    nsplit0 = c['nsplit0']

    if is2n==0:
        nsplit=nsplit0
    else:
        npair_tot0 = get_npair_by_noise(c['s2n_vals'][0], c['desired_err'],c['run'])
        nsplit = int( ceil( nsplit0*float(npair_tot)/npair_tot0 ))
    
    npair_per = int(ceil(npair_tot/float(nsplit)))

    return npair_per, nsplit

def get_npair_nsplit(c, is2n):
    """
    Get number of pairs per split and number of splits

    For equal_time, we take number per split from is2n==0
    """
    if 'desired_err' in c:
        return get_npair_nsplit_by_noise(c, is2n)
    else:
        s2n = c['s2n_vals'][is2n]

        npair = get_s2n_nrepeat(s2n, fac=c['s2n_fac'])
        if npair < c['min_npair']:
            npair = c['min_npair']
        nsplit=c['nsplit']

        return npair,nsplit

