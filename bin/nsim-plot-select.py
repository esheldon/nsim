from __future__ import print_function

import numpy
import fitsio
import biggles
import nsim
import esutil as eu

from argparse import ArgumentParser

parser=ArgumentParser()
parser.add_argument('epsname',help="file in which to write the plot")
parser.add_argument('--stars',action='store_true', help="include stars")
parser.add_argument('--show',action='store_true', help="show the plot")
class DummyArgs(object):
    def __init__(self, runs):
        self.runs=runs

def get_fname(s2n_min, args):
     root='/u/ki/esheldon/lensing/shapesim/run-bdj02mcal02/fit-m-c/run-bdj02mcal02-means-bdj02mcal03-bdj02mcal04-bdj03mcal01-bdj03mcal02-bdj03mcal05-bdj03mcal06-bdj03mcal07-bdj03mcal08-bdj03mcal09-bdj04mcal01-bdj04mcal02-bdj04mcal03'

     if args.stars:
         root += '-bdj03stars-mcal01-bdj03stars-mcal02-bdj03stars-mcal03'

     fname = '%s-select-s2n->-%d.fits' % (root, s2n_min)

     return fname

def read_means(s2n_min, args):
    fname = get_fname(s2n_min, args)

    print('reading:',fname)
    means = fitsio.read(fname)
    means_nocorr = fitsio.read(fname, ext='nocorr')
    return means, means_nocorr

def doplot(s2n_min, data, args):


    """
    data['m'] /= 1.0e-3
    data['merr'] /= 1.0e-3
    data['m_nocorr'] /= 1.0e-3
    data['merr_nocorr'] /= 1.0e-3

    data['c'] /= 1.0e-4
    data['cerr'] /= 1.0e-4
    data['c_nocorr'] /= 1.0e-4
    data['cerr_nocorr'] /= 1.0e-4
    """


    #tab=biggles.Table(2, 1)


    mplt = biggles.FramedPlot()
    mplt.xlabel=r'$(S/N)_{min}$'
    #mplt.ylabel=r'$m [10^{-3}]$'
    mplt.ylabel=r'$m$'
    #mplt.yrange=[-20,12]
    mplt.yrange=[-0.02, 0.02]
    #mplt.x1.draw_ticklabels=False
    mplt.aspect_ratio=1.0/1.618


    mplt.add(
        biggles.FillBetween(
            [5,20], [1e-3,1e-3], 
            [5,20], [-1e-3,-1e-3],
            color='grey90')
    )
    mplt.add( biggles.Curve([5,20],[0,0]))

    
    """
    cplt = biggles.FramedPlot()
    cplt.xlabel=r'$(S/N)_{min}$'
    cplt.ylabel=r'$c [10^{-4}]$'
    cplt.add( biggles.Curve([5,20],[0,0]))
    """

    ccolor='steelblue'
    ncolor='firebrick2'

    psize=2.5


    csyms=[
        {'color':ccolor,'type':'filled circle','size':psize},
        {'color':'black','type':'circle','size':psize},
    ]
    nsyms=[
        {'color':ncolor,'type':'filled diamond','size':psize},
        {'color':'black','type':'diamond','size':psize},
    ]

    clines=[
        {'color':ccolor},
    ]
    nlines=[
        {'color':ncolor},
    ]

    types={
        'mtypes': [
            {
                'name':'m',
                'err':'merr',
                'syms':csyms,
                'lines':clines,
                'label':'corrected for selection',
            }
            ,
            {
                'name':'m_nocorr',
                'err':'merr_nocorr',
                'syms':nsyms,
                'lines':nlines,
                'label':'not corrected',
            }
        ],
    }

    mlist=[]
    for mtype in types['mtypes']:
        for line in mtype['lines']:
            mcurve=biggles.Curve(
                s2n_min,
                data[ mtype['name'] ],
                color=line['color'],
            )
            mplt.add(mcurve)

        for isym,sym in enumerate(mtype['syms']):
            if isym == 0:
                merr=biggles.SymmetricErrorBarsY(
                    s2n_min,
                    data[ mtype['name'] ],
                    data[ mtype['err'] ],
                    color=sym['color'],
                )
                mplt.add(merr)

            mpts=biggles.Points(
                s2n_min,
                data[ mtype['name'] ],
                type=sym['type'],
                size=sym['size'],
                color=sym['color'],
            )
            mplt.add(mpts)

            if isym == 0:
                mpts.label = mtype['label']
                mlist.append( mpts )

    mkey = biggles.PlotKey(0.9, 0.9, mlist, halign='right')
    mplt.add( mkey)

    """
    for ctype in types['ctypes']:
        ccurve=biggles.Curve(
            s2n_min,
            data[ ctype['name'] ],
            color=ctype['color'],
        )
        cpts=biggles.Points(
            s2n_min,
            data[ ctype['name'] ],
            type='filled circle',
            color=ctype['color'],
        )
        cerr=biggles.SymmetricErrorBarsY(
            s2n_min,
            data[ ctype['name'] ],
            data[ ctype['err'] ],
            color=ctype['color'],
        )
        cpts.label = ctype['label']

        cplt.add(ccurve, cpts, cerr)

    tab[0,0] = mplt
    tab[1,0] = cplt
    """

    print("writing:",args.epsname)
    mplt.write_eps(args.epsname)
    if args.show:
        mplt.show()



def main():

    biggles.configure('default','fontsize_min',2.5)

    args=parser.parse_args()

    stars_runs=[
        "run-bdj02mcal02",
        "run-bdj02mcal03",
        "run-bdj02mcal04",
        "run-bdj03mcal01",
        "run-bdj03mcal02",
        "run-bdj03mcal05",
        "run-bdj03mcal06",
        "run-bdj03mcal07",
        "run-bdj03mcal08",
        "run-bdj03mcal09",
        "run-bdj04mcal01",
        "run-bdj04mcal02",
        "run-bdj04mcal03",
    ]
    
    if args.stars:
        runs += [
            "run-bdj03stars-mcal01",
            "run-bdj03stars-mcal02",
            "run-bdj03stars-mcal03",
        ]

    s2n_min = numpy.array( [7,10,13,16,19] )

    fitlist=[]
    for smin in s2n_min:
        means, means_nocorr = read_means(smin, args)

        fits=nsim.averaging_new.get_m_c_oneshear(means)
        fits_nocorr=nsim.averaging_new.get_m_c_oneshear(means_nocorr)

        dt=[
            ('m_nocorr','f8'),
            ('merr_nocorr','f8'),
            ('c_nocorr','f8'),
            ('cerr_nocorr','f8'),
        ]
        allfits=eu.numpy_util.add_fields(fits, dt)
        allfits['m_nocorr'] = fits_nocorr['m']
        allfits['merr_nocorr'] = fits_nocorr['merr']
        allfits['c_nocorr'] = fits_nocorr['c']
        allfits['cerr_nocorr'] = fits_nocorr['cerr']

        fitlist.append(allfits)

    
    data=eu.numpy_util.combine_arrlist(fitlist)

    doplot(s2n_min, data, args)

main()
