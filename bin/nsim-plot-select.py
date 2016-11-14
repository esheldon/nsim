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
     #root='/u/ki/esheldon/lensing/shapesim/run-bdj02mcal02/fit-m-c/run-bdj02mcal02-means-bdj02mcal03-bdj02mcal04-bdj03mcal01-bdj03mcal02-bdj03mcal05-bdj03mcal06-bdj03mcal07-bdj03mcal08-bdj03mcal09-bdj04mcal01-bdj04mcal02-bdj04mcal03'

     if args.stars:
         raise RuntimeError("implement stars runs")
     else:
        fname='/nfs/slac/des/fs1/g/sims/esheldon/lensing/shapesim/runs-bdj03-001/fit-m-c/runs-bdj03-001-means-preselect-select-s2n->-%d.fits' % s2n_min

     #fname = '%s-select-s2n->-%d.fits' % (root, s2n_min)

     return fname

def read_means(s2n_min, args):
    fname = get_fname(s2n_min, args)

    print('reading:',fname)
    means = fitsio.read(fname)
    means_nocorr = fitsio.read(fname, ext='nocorr')
    return means, means_nocorr

def doplot_with_uncorrected(s2n_min, data, args):



    mplt = biggles.FramedPlot()
    mplt.xlabel=r'$(S/N)_{min}$'

    mplt.ylabel=r'$m$'

    mplt.yrange=[-0.02, 0.02]

    mplt.aspect_ratio=1.0/1.618


    mplt.add(
        biggles.FillBetween(
            [5,20], [1e-3,1e-3], 
            [5,20], [-1e-3,-1e-3],
            color='grey90')
    )
    mplt.add( biggles.Curve([5,20],[0,0]))
    

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

    epsname=args.epsname.replace('.eps','-with-nocorr.eps')
    print("writing:",epsname)
    mplt.write_eps(epsname)
    if args.show:
        mplt.show()


def doplot(s2n_min, data_in, args):

    data=data_in.copy()

    data['m'] /= 1.0e-3
    data['merr'] /= 1.0e-3

    data['c'] /= 1.0e-4
    data['cerr'] /= 1.0e-4

    tab=biggles.Table(2, 1)


    mplt = biggles.FramedPlot()
    #mplt.xlabel=r'$(S/N)_{min}$'
    mplt.ylabel=r'$m [10^{-3}]$'
    mplt.yrange=[-2,2]
    mplt.x1.draw_ticklabels=False
    #mplt.aspect_ratio=1.0/1.618


    cmin=5
    cmax=20
    mplt.add(
        biggles.FillBetween(
            [cmin,cmax], [1,1], 
            [cmin,cmax], [-1,-1],
            color='grey90')
    )
    mplt.add( biggles.Curve([cmin,cmax],[0,0]))

    
    cplt = biggles.FramedPlot()
    cplt.xlabel=r'$(S/N)_{min}$'
    cplt.ylabel=r'$c [10^{-4}]$'
    cplt.yrange=[-2,2]
    #cplt.aspect_ratio=1.0/1.618

    cplt.add(
        biggles.FillBetween(
            [cmin,cmax], [1,1], 
            [cmin,cmax], [-1,-1],
            color='grey90')
    )

    cplt.add( biggles.Curve([cmin,cmax],[0,0]))

    color='steelblue'
    sym='filled circle'
    psize=2.5


    csyms=[
        {'color':color,'type':'filled circle','size':psize},
        {'color':'black','type':'circle','size':psize},
    ]

    clines=[
        {'color':color},
    ]

    types={
        'mtypes':
            {'name':'m',
             'err':'merr',
             'syms':csyms,
             'lines':clines,
             'plt':mplt,
            },
        'ctypes':
            {'name':'c',
             'err':'cerr',
             'syms':csyms,
             'lines':clines,
             'plt':cplt,
            },
    }

    for type,types in types.iteritems():
        plt=types['plt']

        for line in types['lines']:
            curve=biggles.Curve(
                s2n_min,
                data[ types['name'] ],
                color=line['color'],
            )
            plt.add(curve)

        for isym,sym in enumerate(types['syms']):
            if isym == 0:
                err=biggles.SymmetricErrorBarsY(
                    s2n_min,
                    data[ types['name'] ],
                    data[ types['err'] ],
                    color=sym['color'],
                )
                plt.add(err)

            pts=biggles.Points(
                s2n_min,
                data[ types['name'] ],
                type=sym['type'],
                size=sym['size'],
                color=sym['color'],
            )
            plt.add(pts)


    tab[0,0] = mplt
    tab[1,0] = cplt

    print("writing:",args.epsname)
    tab.write_eps(args.epsname)
    if args.show:
        tab.show()


def doprint_table(s2n_min, data):

    mfac=1.0/1.0e-3
    cfac=1.0/1.0e-5

    head=r"""
\begin{table*}
    \centering
    \caption{\Mcal\ results for the \bdsim\ simulation with various
        cuts on \snr.   Results are shown with and without corrections
        for selection effects.
    \label{tab:results_sel}}
    \begin{tabular}{ |l| c|c|  c|c|}
        \hline
        & \multicolumn{2}{c}{Uncorrected for Selection}                      & \multicolumn{2}{c}{Corrected for Selection} \\
        Selection                   & $m$             & $c$            & $m$               & $c$  \\
                                    & $[10^{-3}]$     & $[10^{-5}]$    & $[10^{-3}]$       & $[10^{-5}]$ \\
        \hline"""
        
    body=[]
    fmt='%+.2f'
    efmt='%.2f'
    formats=['$%(fmt)s \\pm %(efmt)s$' % dict(fmt=fmt,efmt=efmt)]*4
    formats = ' & '.join(formats)
    pattern=r"$\mbox{\snr} > %d $ & " + formats + r" \\"
    for i,s2n in enumerate(s2n_min):
        m=data['m'][i]*mfac
        merr=data['merr'][i]*mfac
        m_nocorr=data['m_nocorr'][i]*mfac
        merr_nocorr=data['merr_nocorr'][i]*mfac

        c=data['c'][i]*cfac
        cerr=data['cerr'][i]*cfac
        c_nocorr=data['c_nocorr'][i]*cfac
        cerr_nocorr=data['cerr_nocorr'][i]*cfac

        line=pattern % (
            s2n,
            m_nocorr,merr_nocorr,
            c_nocorr,cerr_nocorr,
            m,merr,
            c,cerr,
        )
        body.append(line)


    tail=r"""    \end{tabular}
\end{table*}
"""
    
    table = [head] + body + [tail]

    table='\n'.join(table)
    
    print(table)

def main():

    biggles.configure('default','fontsize_min',2.5)

    args=parser.parse_args()

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

    doplot_with_uncorrected(s2n_min, data, args)
    doplot(s2n_min, data, args)

    doprint_table(s2n_min, data)

main()
