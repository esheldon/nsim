from sys import stderr
import os
from os.path import join as path_join

import yaml
import numpy
from numpy import array


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
            f += '-%05i' % itrial
    f += '.%s' % ext
    return path_join(dir, f)

def read_output(run, is2n, fs=None, ext='fits'):
    """
    Read the collated file with all trials
    """
    import fitsio
    fname=get_output_url(run, 0, is2n, fs=fs, ext=ext)
    return fitsio.read(fname)


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


# from eg01r04 with bugged code
s2n_ref_eg_sr2=array([ 10, 15, 23, 35, 53, 81, 123, 187, 285, 433, 658, 1000 ],dtype='f8')
err_ref_eg_sr2=array([5.86094527e-05,   4.24659823e-05,   4.15026596e-05,
                      4.10366038e-05,   4.17969006e-05,   4.19735634e-05,
                      4.19306727e-05,   4.19540700e-05,   4.14286130e-05,
                      3.82807253e-05,   2.54193529e-05,   1.65443563e-05])

npair_ref_eg_sr2=array([2422000, 2422000, 1175708,  539933,  229964,   98212,   43000,
                        18684,    8160,    4150,    4150,    4150])

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
            npairii = numpy.interp([s2n], s2n_ref_eg_sr2, npair_ref_eg_sr2)
            errii = numpy.interp([s2n], s2n_ref_eg_sr2, err_ref_eg_sr2)
        elif sigratio >= 1.41:
            npairii = numpy.interp([s2n], s2n_ref_eg_sr14, npair_ref_eg_sr14)
            errii = numpy.interp([s2n], s2n_ref_eg_sr14, err_ref_eg_sr14)
        else:
            npairii = numpy.interp([s2n], s2n_ref_eg_sr1, npair_ref_eg_sr1)
            errii = numpy.interp([s2n], s2n_ref_eg_sr1, err_ref_eg_sr1)
    elif 'bdfg' in run:
        # for bdf make sure you set sec_per_pair
        npairii = numpy.interp([s2n], s2n_ref_bdfg, npair_ref_bdfg)
        errii = numpy.interp([s2n], s2n_ref_bdfg, err_ref_bdfg)
    elif '-dg' in run:
        if sigratio >= 1.99:
            npairii = numpy.interp([s2n], s2n_ref_dg_sr2, npair_ref_dg_sr2)
            errii = numpy.interp([s2n], s2n_ref_dg_sr2, err_ref_dg_sr2)
        elif sigratio >= 1.41:
            npairii = numpy.interp([s2n], s2n_ref_dg_sr14, npair_ref_dg_sr14)
            errii = numpy.interp([s2n], s2n_ref_dg_sr14, err_ref_dg_sr14)
        else:
            npairii = numpy.interp([s2n], s2n_ref_dg_sr1, npair_ref_dg_sr1)
            errii = numpy.interp([s2n], s2n_ref_dg_sr1, err_ref_dg_sr1)
    elif '-gg' in run:
        if sigratio >= 1.99:
            npairii = numpy.interp([s2n], s2n_ref_gg_sr2, npair_ref_gg_sr2)
            errii = numpy.interp([s2n], s2n_ref_gg_sr2, err_ref_gg_sr2)
        elif sigratio >= 1.41:
            npairii = numpy.interp([s2n], s2n_ref_gg_sr14, npair_ref_gg_sr14)
            errii = numpy.interp([s2n], s2n_ref_gg_sr14, err_ref_gg_sr14)
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

    if 'obj_T_mean' not in c['simc']:
        # for now; this might be conservative?
        sigratio=2.0
    else:
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

