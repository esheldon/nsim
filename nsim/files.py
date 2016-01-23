from __future__ import print_function
import os
from os.path import join as path_join

import yaml
import numpy
from numpy import array, sqrt, log10


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


#
# wq
#

def get_wq_dir(run):
    dir=get_run_dir(run)
    dir=path_join(dir, 'wq')
    return dir

def get_wq_job_url(run, filenum, missing=False):
    d=get_wq_dir(run)

    end = '-%06d' % filenum
    if missing:
        end='%s-missing' % end

    fname='{run}{end}.yaml'.format(run=run,end=end)
    return path_join(d,fname)

def get_wq_combine_psample_job_url(run, is2n, itrial):
    d=get_wq_dir(run)

    fname='%s-%05d-%05d-comb-psample.yaml' % (run,is2n,itrial)
    return path_join(d,fname)

def get_wq_master_url(run):
    d=get_wq_dir(run)
    return path_join(d,'%s.sh' % run)

#
# lsf
#

def get_lsf_dir(run):
    dir=get_run_dir(run)
    dir=path_join(dir, 'lsf')
    return dir

def get_lsf_job_url(run, filenum, missing=False):
    d=get_lsf_dir(run)

    end = '-%06d' % filenum
    if missing:
        end='%s-missing' % end

    fname='{run}{end}.lsf'.format(run=run,end=end)
    return path_join(d,fname)

def get_lsf_combine_psample_job_url(run, is2n, itrial):
    d=get_lsf_dir(run)

    fname='%s-%05d-%05d-comb-psample.lsf' % (run,is2n,itrial)
    return path_join(d,fname)


def get_lsf_master_url(run):
    d=get_lsf_dir(run)
    return path_join(d,'%s.sh' % run)


def get_slr_dir(run):
    dir=get_run_dir(run)
    dir=path_join(dir, 'slr')
    return dir

def get_slr_job_url(run):
    d=get_slr_dir(run)

    fname='%s.slr' % run
    return path_join(d,fname)

#def get_slr_master_url(run):
#    d=get_slr_dir(run)
#    return path_join(d,'%s.sh' % run)


def get_output_dir(run, sub=None, fs=None):
    dir=get_run_dir(run, fs=fs)
    dir=path_join(dir, 'outputs')
    if sub:
        dir = path_join(dir, sub)
    return dir

def get_output_url(run, is2, is2n, itrial=None, fs=None, ext='fits'):
    """
    is2 and is2n are the index in the list of s2 and s2n vals for a given run.
    """
    sub=None
    if itrial is not None:
        sub='bytrial'
    dir=get_output_dir(run, sub=sub, fs=fs)
    f='%s-%03i-%03i' % (run,is2,is2n)
    if itrial is not None:
        if itrial == '*':
            f += '-*'
        else:
            f += '-%06d' % itrial
    f += '.%s' % ext
    return path_join(dir, f)

def get_noise_dir(run):
    dir=get_run_dir(run)
    dir=path_join(dir, 'noise')
    return dir

def get_noise_file(run):
    dir=get_noise_dir(run)
    fname='%s-noise-im.fits'
    fname=os.path.join(dir, fname)
    return fname

def get_plot_dir(run):
    dir=get_run_dir(run)
    dir=path_join(dir, 'plots')
    return dir

def get_plot_url(run, extra, ext='eps'):
    dir=get_plot_dir(run)
    f='%s-%s.%s' % (run,extra,ext)
    return path_join(dir, f)

def get_means_dir(run):
    dir=get_run_dir(run)
    dir=path_join(dir, 'fit-m-c')
    return dir

def get_means_url(run, extra=None):
    dir=get_means_dir(run)
    if extra is not None:
        extra = '-'+extra
    else:
        extra=''

    f='%s-means%s.fits' % (run,extra)
    return path_join(dir, f)

def get_psample_summed_url(run, is2n, itrial=None, fs=None, ext='fits'):
    """
    is2 and is2n are the index in the list of s2 and s2n vals for a given run.
    """

    sub=None
    if itrial is not None:
        sub='bytrial'

    dir=get_output_dir(run, sub=sub, fs=fs)

    f='%s-%03i-%03i' % (run,0,is2n)

    if itrial is not None:
        if itrial == '*':
            f += '-*'
        else:
            f += '-%05i' % itrial

    f += '-summed.%s' % ext
    return path_join(dir, f)


def read_output(run, is2n, fs=None, itrial=None, ext='fits', **kw):
    """
    Read the collated file with all trials
    """
    import fitsio
    fname=get_output_url(run, 0, is2n, itrial=itrial, fs=fs, ext=ext)
    print("reading collated file:",fname)
    return fitsio.read(fname, **kw)

def get_fitprior_url(run, is2n, itrial=None, fs=None, extra=None, ext='fits'):
    """
    we fit prior from high s/n sample
    """
    url=get_output_url(run, 0, is2n, itrial=itrial, fs=fs, ext=ext)

    end=['fitprior']
    if extra is not None:
        end.append(extra)

    end = '-'.join(end)

    url = url.replace('.%s' % ext,'-%s.%s' % (end,ext))
    return url

def get_extra_url(bname):
    """
    return the url for the file named {simdir}/bname
    """
    d=get_simdir()
    f=os.path.join(d,'extra_files',bname)

    return f

def get_averaged_url(run, is2n=None, fs=None, ext='fits'):
    """
    send is2n= for the split one
    """

    dir=get_output_dir(run, fs=fs)
    if is2n is not None:
        f='%s-%03i-%03i-avg' % (run,0,is2n)
    else:
        f='%s-%03i-avg' % (run,0)
    f = '%s.%s' % (f,ext)
    return path_join(dir, f)

def read_averaged(run, is2n=None, fs=None, ext='fits'):
    """
    Read the file with all averaged and summed quantities for
    each s/n bin
    """
    import fitsio
    fname=get_averaged_url(run, is2n=is2n, fs=fs, ext=ext)
    return fitsio.read(fname)


s2n_ref_bdfg=[ 10, 15, 23, 35, 53, 81, 123, 187, 285, 433, 658, 1000 ]
npair_ref_bdfg=[1240000, 1240000,  711574,  363878,  164300,
               71517,  31589,  13733,   6045,
               3038,   3038,   3038]
err_ref_bdfg=[ 8.64686523781e-05,6.36021322793e-05,5.78428340247e-05,5.44774877257e-05,5.41502412604e-05,5.39701794399e-05,5.35799931151e-05,5.38147319202e-05,5.35173339764e-05,4.97357492734e-05,3.24657802612e-05,2.14956557248e-05]


# for flux lim sigratio > 0.85
s2n_ref_eg_fluxlim85=array([ 5.0, 10.7, 23., 48., 103., 220., 469., 1000.0 ],
                           dtype='f8')
err_ref_eg_fluxlim85=array([3.77692058e-05,   4.30906570e-05,   4.15750598e-05,
                            4.13269876e-05,   4.70262497e-05,   5.09335013e-05,
                            6.04725966e-05,   4.16510351e-05],
                           dtype='f8')
npair_ref_eg_fluxlim85=array([9700000, 3080138,  811114,  219026,
                              48500,   10615,    2184, 1002],dtype='f8')

# from eg01r04 with bugged code
# 5 is fake right now
s2n_ref_eg_sr2=array([5, 10, 15, 23, 35, 53, 81, 123, 187, 285, 433, 658, 1000 ],dtype='f8')
err_ref_eg_sr2=array([5.0e-05,5.86094527e-05,   4.24659823e-05,   4.15026596e-05,
                      4.10366038e-05,   4.17969006e-05,   4.19735634e-05,
                      4.19306727e-05,   4.19540700e-05,   4.14286130e-05,
                      3.82807253e-05,   2.54193529e-05,   1.65443563e-05])
npair_ref_eg_sr2=array([9688000, 2422000, 2422000, 1175708,  539933,  229964,   98212,   43000,
                        18684,    8160,    4150,    4150,    4150])

#s2n_ref_eg_sr2=array([ 10, 15, 23, 35, 53, 81, 123, 187, 285, 433, 658, 1000 ],dtype='f8')
#err_ref_eg_sr2=array([5.86094527e-05,   4.24659823e-05,   4.15026596e-05,
#                      4.10366038e-05,   4.17969006e-05,   4.19735634e-05,
#                      4.19306727e-05,   4.19540700e-05,   4.14286130e-05,
#                      3.82807253e-05,   2.54193529e-05,   1.65443563e-05])
#npair_ref_eg_sr2=array([2422000, 2422000, 1175708,  539933,  229964,   98212,   43000,
#                        18684,    8160,    4150,    4150,    4150])

s2n_ref_eg_sr14 = array([ 10, 15, 23, 35, 53, 81, 123, 187, 285, 433, 658, 1000 ],dtype='f8')
err_ref_eg_sr14 = array([  7.03558581e-05,   5.31578898e-05,   5.33730462e-05,
                         5.33931400e-05,   5.45184327e-05,   5.49625225e-05,
                         5.48505289e-05,   5.49522236e-05,   5.43305311e-05,
                         5.05090984e-05,   3.33051218e-05,   2.15534084e-05])


npair_ref_eg_sr14= array([2420000, 2420000, 1174789,  539418,  230989,   98736,   43197,
                          18600,    8160,    4165,    4165,    4165])


#s2n_ref_eg_sr1 = array([ 10, 15, 23, 35, 53, 81, 123, 187, 285, 433, 658, 1000 ],dtype='f8')
#err_ref_eg_sr1 = array([8.39715840e-05,   6.76302528e-05,   7.12956923e-05,
#                        7.33158683e-05,   7.59915389e-05,   7.68021383e-05,
#                        7.65095703e-05,   7.72849329e-05,   7.62656351e-05,
#                        7.01269993e-05,   4.64553805e-05,   3.04000144e-05])

#npair_ref_eg_sr1=array([2412000, 2412000, 1171026,  537675,  230346,   98490,   43014,
#                        18600,    8159,    4158,    4158,    4158])

# s/n 5 is fake right now
s2n_ref_eg_sr1 = array([ 5, 10, 15, 23, 35, 53, 81, 123, 187, 285, 433, 658, 1000 ],dtype='f8')
err_ref_eg_sr1 = array([9.0e-5, 8.39715840e-05,   6.76302528e-05,   7.12956923e-05,
                        7.33158683e-05,   7.59915389e-05,   7.68021383e-05,
                        7.65095703e-05,   7.72849329e-05,   7.62656351e-05,
                        7.01269993e-05,   4.64553805e-05,   3.04000144e-05])

npair_ref_eg_sr1=array([8412000,2412000, 2412000, 1171026,  537675,  230346,   98490,   43014,
                        18600,    8159,    4158,    4158,    4158])


s2n_ref_eg_sr084 = array([ 5, 11, 22, 47, 100 ], dtype='f8')

err_ref_eg_sr084 =array([ 0.00038262,  0.00036024,  0.00027844,  0.00024308,  0.0002513 ])
npair_ref_eg_sr084=array([1703550,  392850,  158400,   41550,    9600])


# from dg01r01 jackknifed
s2n_ref_dg_sr2=array([ 10, 15, 23, 35, 53, 81, 123, 187, 285, 433, 658, 1000 ],dtype='f8')

err_ref_dg_sr2= array([5.87276281e-05,   4.88444115e-05,   4.45253939e-05,
                       4.24857060e-05,   4.20594285e-05,   4.22970136e-05,
                       4.17526066e-05,   4.15213525e-05,   4.30400294e-05,
                       4.23330372e-05,   4.37362379e-05,   4.22950929e-05])

npair_ref_dg_sr2= array([14313377,  9448873,  5179365,  2546812,  1180505,   518205,
                         226594,    98256,    42364,    18585,     8010,     3460])


# from dg04r02 jackknifed
s2n_ref_dg_sr14=array([ 10, 15, 23, 35, 53, 81, 123, 187, 285, 433, 658, 1000 ],dtype='f8')

err_ref_dg_sr14= array([9.36238578e-05,   7.45277261e-05,   6.50421446e-05,
                        6.06075219e-05,   6.02431965e-05,   5.95502050e-05,
                        5.96811389e-05,   6.22620554e-05,   5.89997635e-05,
                        6.06547374e-05,   5.84964467e-05,   6.02839043e-05])

npair_ref_dg_sr14= array([8200000, 5849224, 3446460, 1776776,  845748,  377036,  166132,
                          72488,   31652,   13776,    5920,    2512])



# with jackknife from run-dg05r01
s2n_ref_dg_sr1=array([ 10, 15, 23, 35, 53, 81, 123, 187, 285, 433, 658, 1000 ],dtype='f8')

err_ref_dg_sr1= array( [0.00010934769005102428, 8.34675545843533e-05, 6.8885858378181809e-05,
                        6.2408171398520217e-05, 6.0111501604497353e-05, 5.9265828120324162e-05,
                        5.94227339673447e-05, 5.9275576860331574e-05, 5.927755979925429e-05,
                        6.0059247335501772e-05, 6.0083958777380203e-05, 
                        5.9532082538788167e-05])
    
npair_ref_dg_sr1= array([9279768, 7184112, 4646496, 2575200, 1278552, 582784,
                         259144, 113216, 48720, 21021, 9200, 4032])


# old one without jackknife
#s2n_ref_dg_sr1=array([ 10, 15, 23, 35, 53, 81, 123, 187, 285, 433, 658, 1000 ],dtype='f8')

#err_ref_dg_sr1= array([6.15830895e-05,   5.41834995e-05,   5.75238159e-05,
#                       5.98859154e-05,   6.27954858e-05,   6.42454388e-05,
#                       6.44668215e-05,   6.45949509e-05,   6.39895879e-05,
#                       5.94180555e-05,   3.92390091e-05,   2.59813120e-05])

#npair_ref_dg_sr1= array([6100000, 6100000, 3500424, 1789984,  808250,  351848,  155428,
#                         67588,   29646,   14884,   14884,   14884])



# from gg01r01
s2n_ref_gg_sr2=array([ 10, 15, 23, 35, 53, 81, 123, 187, 285, 433, 658, 1000 ],dtype='f8')
err_ref_gg_sr2=array([4.98680203e-05,   4.98357873e-05,   4.99109235e-05,
                      4.98672236e-05,   4.98538152e-05,   5.00366767e-05,
                      4.97837537e-05,   4.97515484e-05,   4.98803397e-05,
                      4.99418567e-05,   4.95626495e-05,   4.99953840e-05])
npair_ref_gg_sr2 = array([2200000, 1084160,  482895,  213525,   93732,   40296,   17440,
                          7560,    3312,    1421,     612,     266])

s2n_ref_gg_sr14 = array([ 10, 15, 23, 35, 53, 81, 123, 187, 285, 433, 658, 1000 ],dtype='f8')
err_ref_gg_sr14 = array([5.61174116e-05,   5.80381439e-05,   5.90672346e-05,
                         5.94514283e-05,   5.96208412e-05,   6.00325336e-05,
                         5.97678990e-05,   6.01893373e-05,   5.92819589e-05,
                         6.06123456e-05,   5.93288780e-05,   6.02809294e-05])
npair_ref_gg_sr14 =array([2720000, 1340416,  599624,  265064,  116416,   49680,   21600,
                          9315,    4092,    1755,     756,     330])


# from gg04r01
s2n_ref_gg_sr1=array([ 10, 15, 23, 35, 53, 81, 123, 187, 285, 433, 658, 1000 ],dtype='f8')
err_ref_gg_sr1=array([3.52809767e-05,   3.53038557e-05,   3.53300171e-05,
                      3.53451981e-05,   3.52917768e-05,   3.52483577e-05,
                      3.52699247e-05,   3.50577662e-05,   3.49212158e-05,
                      3.45652208e-05,   3.32411166e-05,   3.13146698e-05])

npair_ref_gg_sr1 = array([11050000,  6586463,  3311243,  1536480,   692560,   299420,
                          130240,    57200,    24860,    11016,     5160,     2484])


def get_npair_by_noise(s2n, desired_err, run, sigratio):
    """
    given the desired final error, determine the required number of pairs
    """

    if '-eg' in run:
        if sigratio >= 1.99:
            s2n_ref = s2n_ref_eg_sr2
            npair_ref = npair_ref_eg_sr2
            err_ref = err_ref_eg_sr2
        elif sigratio >= 1.41:
            s2n_ref = s2n_ref_eg_sr14
            npair_ref = npair_ref_eg_sr14
            err_ref = err_ref_eg_sr14
        elif sigratio >= 1.0:
            s2n_ref = s2n_ref_eg_sr1
            npair_ref = npair_ref_eg_sr1
            err_ref = err_ref_eg_sr1
        else:
            s2n_ref = s2n_ref_eg_sr084
            npair_ref = npair_ref_eg_sr084
            err_ref = err_ref_eg_sr084

    elif 'bdfg' in run:
        # for bdf make sure you set sec_per_pair
        s2n_ref = s2n_ref_bdfg
        npair_ref = npair_ref_bdfg
        err_ref = err_ref_eg_bdfg
    elif '-dg' in run:
        if sigratio >= 1.99:
            s2n_ref = s2n_ref_dg_sr2
            npair_ref = npair_ref_dg_sr2
            err_ref = err_ref_dg_sr2
        elif sigratio >= 1.41:
            s2n_ref = s2n_ref_dg_sr14
            npair_ref = npair_ref_dg_sr14
            err_ref = err_ref_dg_sr14
        else:
            s2n_ref = s2n_ref_dg_sr1
            npair_ref = npair_ref_dg_sr1
            err_ref = err_ref_dg_sr1

    elif '-gg' in run:
        if sigratio >= 1.99:
            s2n_ref = s2n_ref_gg_sr2
            npair_ref = npair_ref_gg_sr2
            err_ref = err_ref_gg_sr2
        elif sigratio >= 1.41:
            s2n_ref = s2n_ref_gg_sr14
            npair_ref = npair_ref_gg_sr14
            err_ref = err_ref_gg_sr14
        else:
            s2n_ref = s2n_ref_gg_sr1
            npair_ref = npair_ref_gg_sr1
            err_ref = err_ref_gg_sr1

    else:
        raise ValueError("support runs of type '%s'" % run)

    log_s2n=log10( s2n_ref )
    log_npair=log10( npair_ref )

    log_npairii = numpy.interp(log10([s2n]), log_s2n, log_npair)
    npairii = 10.0**log_npairii

    errii = numpy.interp([s2n], s2n_ref, err_ref)


    # desired_err = errii*sqrt(npairii/npair)
    # thus npair = npairii*(errii/desired_err)^2
    npair = npairii*(errii/desired_err)**2

    return npair[0]
    
def get_npair_by_noise_fluxlim(s2n, desired_err, run, sigratio_low):
    """
    given the desired final error, determine the required number of pairs
    """

    if '-eg' in run:
        print('DOING FLUXLIM')
        if sigratio_low < 0.85:

            s2n_ref   = s2n_ref_eg_fluxlim85 
            npair_ref = npair_ref_eg_fluxlim85 
            err_ref   = err_ref_eg_fluxlim85 

    else:
        raise ValueError("support '%s'" % run)


    log_s2n=log10( s2n_ref )
    log_npair=log10( npair_ref )

    log_npairii = numpy.interp(log10([s2n]), log_s2n, log_npair)
    npairii = 10.0**log_npairii

    errii = numpy.interp([s2n], s2n_ref, err_ref)

    npair = npairii*(errii/desired_err)**2
    print('NPAIR for s/n=%s is %s' % (s2n,npair) )
    return npair

def get_npair_nsplit(c, is2n, npair_min=None):
    """
    Get the npair/nsplit given the input config and s2n index
    """
    from math import ceil

    if 'npair' in c:
        npair_tot=c['npair'][is2n]

        # this is per
        tmsec = c['desired_hours']*3600.0
        seconds_per=c['sec_per_pair']
        npair_per = int(round( tmsec/seconds_per ))

        nsplit = int(ceil( npair_tot/float(npair_per) ))

        return npair_per, nsplit



    # this is the requirement from measurement error
    s2n = c['s2n_vals'][is2n]

    if 'joint_TF_dist' in c['simc']:

        #
        # This is the lowest sig ratio now
        #

        sigratio=sqrt(c['simc']['T_bounds'][0]/c['simc']['psf_T'])
        npair_tot = get_npair_by_noise_fluxlim(s2n,
                                               c['desired_err'],c['run'],
                                               sigratio)
    else:
        simc=c['simc']
        if 'obj_T' in simc:
            T=simc['obj_T']
        else:
            T=simc['obj_T_mean']

        sigratio = sqrt( T/simc['psf_T'] )
        npair_tot = get_npair_by_noise(s2n, c['desired_err'],c['run'], sigratio)

    npair_shapenoise = 0
    do_ring=simc.get('do_ring',True)
    if not do_ring:
        print("not doing ring")
        # add in shape noise term for non-ring test
        ngal = (0.22/c['desired_err'])**2
        npair_shapenoise = int(ngal)

    npair_tot += npair_shapenoise

    if npair_min is not None and npair_tot < npair_min:
        npair_tot = npair_min

    # desired time per split
    # = sec_per_pair * npair_per_split

    # convert from hours to seconds
    tmsec = c['desired_hours']*3600.0
    
    npair_per = tmsec/c['sec_per_pair']
    npair_per = int(ceil(npair_per))
    nsplit = int(ceil( npair_tot/float(npair_per) ))

    return npair_per, nsplit

def get_gal_nsplit(c):
    """
    new non ring where we request a specific number of gals
    """
    from math import ceil

    ngal = c['ngal']
    nrand = c.get('nrand',1)

    # multiply by nrand to get total time per galaxy
    tmsec = c['desired_hours']*3600.0

    sec_per = c['sec_per']*nrand

    ngal_per = int(round( tmsec/sec_per ) )

    nsplit = int(ceil( ngal/float(ngal_per) ))

    time_hours = ngal_per*sec_per/3600.0

    print("ngal requested:",ngal,"nrand:",nrand)
    print('seconds per image:',c['sec_per'],"sec per with rand:",sec_per)
    print('nsplit:',nsplit,'ngal per:',ngal_per,'time (hours):',time_hours)


    return ngal_per, nsplit, time_hours
