from __future__ import print_function

class DiscreteSampler(object):
    def __init__(self, sample_dict):
        """
        sample_dict is 
        {key1: frequency1, key2:frequency2, ...}

        Where the frequencies are arbitrary, can be in [0,1] or [0,100] etc.
        """
        self.sample_dict=sample_dict
        self.wsum = float( sum(self.sample_dict.values()) )

    def __call__(self):
        import random 
        number = random.random() * self.wsum
        for k,v in self.sample_dict.iteritems():
            if number < v:
                break
            number -= v
        return k

def test(nsample=10000, ntypes=2):
    # the following values can be any non-negative numbers, no need of sum=100
    if ntypes==2:
        weights = {'set1': 90.0,
                   'set2': 10.0}
    elif ntypes==3:
        weights = {'set1': 80,
                   'set2': 15,
                   'set3': 5}
    else:
        raise ValueError("ntypes=2 or 3")

    ds=DiscreteSampler(weights)

    samples = [ds() for i in xrange(nsample)]

    
    print("type     truefreq         freq")
    print("------------------------------")
    for key in weights:
        freq = samples.count(key)/float(nsample)
        print("%s %12g %12g" % (key,weights[key]/ds.wsum,freq))

    return ds
