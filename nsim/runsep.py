import logging
import numpy as np
logger = logging.getLogger(__name__)

SEP_FILTER_KERNEL =  np.array([
    [0.004963, 0.021388, 0.051328, 0.068707, 0.051328, 0.021388, 0.004963],
    [0.021388, 0.092163, 0.221178, 0.296069, 0.221178, 0.092163, 0.021388],
    [0.051328, 0.221178, 0.530797, 0.710525, 0.530797, 0.221178, 0.051328],
    [0.068707, 0.296069, 0.710525, 0.951108, 0.710525, 0.296069, 0.068707],
    [0.051328, 0.221178, 0.530797, 0.710525, 0.530797, 0.221178, 0.051328],
    [0.021388, 0.092163, 0.221178, 0.296069, 0.221178, 0.092163, 0.021388],
    [0.004963, 0.021388, 0.051328, 0.068707, 0.051328, 0.021388, 0.004963],
])

SEP_THRESH=0.8
SEP_PARS={
    'deblend_cont':0.00001,
    'deblend_nthresh':64,
    'minarea':4,
    'filter_kernel':SEP_FILTER_KERNEL,
}

def find_objects(obs):
    import sep

    noise=np.sqrt(1.0/obs.weight[0,0])
    objs=sep.extract(
        obs.image,
        SEP_THRESH,
        err=noise,
        **SEP_PARS
    )
    return objs
