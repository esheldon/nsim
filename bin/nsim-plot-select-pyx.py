from __future__ import print_function
from numpy import array
from pyx import *
from pyx.graph import axis

import pyxtools
from pyxtools import colors, linestyles

# xscale scales everything by that factor, such as
# text, ticks etc.
# vscale: used for sizes of visual elements, like arrows, graph symbols, axis ticks, etc.
# wscale: widths such as line widths

unit.set(xscale=1.25, vscale=1.4)

thresh=[0,7,10,13,16,19]
#s2n_means= [8.68, 15.47, 97.44]
s2n_means = [10.3982, 24.0881, 55.41, 127.428]

xrng = [-2, 21]
log_xrng=[7, 200]

symsize=0.1

corr_color=colors('steelblue')
#nocorr_color=colors('SpringGreen4')
nocorr_color=colors('firebrick2')

def doplot(plt,xvals, yvals, yerr, yvals_corr, yerr_corr):

    corr_values=graph.data.values(
        x=list(xvals),
        y=list(yvals_corr),
        dy=list(yerr_corr),
        title='corrected for selection',
    )

    corr_symbol=graph.style.symbol(
        symbol=graph.style.symbol.circle,
        size=symsize,
        symbolattrs=[deco.filled([corr_color])],
        #symbolattrs=[corr_color],
    )
    corr_line = graph.style.line(
        lineattrs=[corr_color],
    )
    corr_errorbars = graph.style.errorbar(errorbarattrs=[corr_color])
    corr_styles=[corr_line,corr_errorbars,corr_symbol]
    plt.plot(corr_values,corr_styles)

    nocorr_values=graph.data.values(
        x=list(xvals),
        y=list(yvals),
        dy=list(yerr),
        title='not corrected',
    )

    nocorr_symbol=graph.style.symbol(
        symbol=graph.style.symbol.triangle,
        size=symsize,
        symbolattrs=[deco.filled([nocorr_color])],
        #symbolattrs=[nocorr_color],
    )
    nocorr_line = graph.style.line(
        lineattrs=[nocorr_color,style.linestyle.dashed],
    )
    nocorr_errorbars = graph.style.errorbar(errorbarattrs=[nocorr_color])

    nocorr_styles=[nocorr_line,nocorr_errorbars,nocorr_symbol]
    plt.plot(nocorr_values,nocorr_styles)


def doplot_full(x, corr, err_corr, nocorr, err_nocorr, xaxis, yaxis, ypos, key):

    c=canvas.canvas()

    g = graph.graphxy(
        width=8,
        ypos=ypos,
        x=xaxis,
        y=yaxis,
        key=key,
    )

    # need to draw the box first
    #x1,y1 = g.pos(xrng[0], -1)
    #x2,y2 = g.pos(xrng[1], 1)
    if isinstance(xaxis,graph.axis.linkedaxis):
        xmin=xaxis.axis.min
        xmax=xaxis.axis.max
    else:
        xmin=xaxis.min
        xmax=xaxis.max

    x1,y1 = g.pos(xmin, -1)
    x2,y2 = g.pos(xmax, 1)
    b = box.rect(x1, y1, x2-x1, y2-y1)
    c.draw(b.path(),[deco.filled([color.gray(0.85)])])

    doplot(
        g,
        x,
        nocorr, err_nocorr,
        corr, err_corr,
    )
    g.plot(graph.data.function("y(x)=0",title=None))

    c.insert(g)
    return c,g

def plot_m(x, corr, err_corr, nocorr, err_nocorr, ymin, ymax, cplt):

    xaxis=graph.axis.linkedaxis(cplt.axes["x"])

    ylabel=r'$m~[10^{-3}]$'
    yaxis=axis.lin(min=ymin,max=ymax,title=ylabel,density=2)

    ypos=cplt.height+0.5
    key=graph.key.key(pos="br")

    c,g = doplot_full(x,corr,err_corr,nocorr,err_nocorr,xaxis,yaxis,ypos,key)
    return c,g

def plot_c(x, corr, err_corr, nocorr, err_nocorr, ymin, ymax, xaxis):

    ylabel=r'$c~[10^{-4}]$'
    yaxis=axis.lin(min=ymin,max=ymax,title=ylabel,density=2)

    ypos=0.0
    key=None

    c,g = doplot_full(x,corr,err_corr,nocorr,err_nocorr,xaxis,yaxis,ypos,key)
    return c,g

def plot_mrng(cplt):

    c=canvas.canvas()

    ymin,ymax=-6,6

    #nocorr=[-0.37, 7.44, 0.47]
    #err_nocorr=[0.55, 0.37, 0.23]

    #corr=[-0.11, -0.41, -0.02]
    #err_corr=[0.55, 0.37, 0.23]

    nocorr=array([ 3.676e-03,
                  4.550e-03,
                  2.095e-04,
                  -1.998e-03])/1.0e-3
    err_nocorr=array([3.817e-04,
                      3.341e-04,
                      3.717e-04,
                      5.190e-04])/1.0e-3

    corr=array([-3.725e-04,
                -4.067e-04,
                1.429e-04,
                2.215e-04])/1.0e-3
    err_corr=array([3.802e-04,
                    3.325e-04,
                    3.717e-04,
                    5.202e-04 ])/1.0e-3

    c,g=plot_m(s2n_means, corr, err_corr, nocorr, err_nocorr, ymin, ymax, cplt)

    return c,g


def plot_mthresh(cplt):

    c=canvas.canvas()

    #ymin,ymax=-5,5
    ymin,ymax=-6,6
    corr=[-0.91, -0.15, -0.10, -0.09, -0.04, -0.23]
    err_corr=[0.19,0.19,0.20,0.20,0.21,0.22]

    nocorr=[-0.91, 2.53, 3.48, 2.61, 1.85, 0.98]
    err_nocorr=[0.19,0.19,0.20,0.20,0.21,0.22]

    c,g=plot_m(thresh, corr, err_corr, nocorr, err_nocorr, ymin, ymax, cplt)

    return c,g


def plot_c2rng():
    #nocorr     = [1.78, 1.16, 0.46]
    #err_nocorr = [0.27, 0.18, 0.12]

    #corr       = [0.27, 0.09, 0.20]
    #err_corr   = [0.27, 0.18, 0.12]

    nocorr=array([1.702e-04,
                  7.042e-05,
                  4.098e-06,
                  8.133e-05])/1.0e-4

    err_nocorr=array([1.899e-05,
                      1.662e-05,
                      1.849e-05,
                      2.582e-05])/1.0e-4

    corr=array([2.279e-05,
                2.110e-05,
                -1.680e-05,
                6.090e-05])/1.0e-4

    err_corr=array([1.891e-05,
                    1.654e-05,
                    1.849e-05,
                    2.588e-05])/1.0e-4

    ymin,ymax=-3,3
    xlabel=r'$\langle S/N \rangle$'
    xaxis=axis.log(min=log_xrng[0],max=log_xrng[1],title=xlabel)

    c,g=plot_c(s2n_means,corr,err_corr,nocorr,err_nocorr,ymin,ymax,xaxis)
    return c,g



def plot_c2thresh():

    nocorr     = [0.17, 0.93, 0.81, 0.58, 0.49, 0.42]
    err_nocorr = [0.09, 0.10, 0.10, 0.10, 0.11, 0.11]

    corr       = [0.17, 0.17, 0.15, 0.13, 0.14, 0.12]
    err_corr   = [0.09, 0.10, 0.10, 0.10, 0.11, 0.11]
    xlabel=r'$S/N$ threshold'
    ymin,ymax=-3,3

    xaxis=axis.lin(min=xrng[0],max=xrng[1],title=xlabel)

    c,g=plot_c(thresh,corr,err_corr,nocorr,err_nocorr,ymin,ymax,xaxis)
    return c,g


def do_all_plots():
    cth_plt, cth_g=plot_c2thresh()
    mth_plt, mth_g=plot_mthresh(cth_g)

    thresh_canvas = canvas.canvas()
    thresh_canvas.insert(cth_plt)
    thresh_canvas.insert(mth_plt)
    thresh_canvas.writetofile('mc-select-bias-thresh.eps')
    thresh_canvas.writetofile('mc-select-bias-thresh.pdf')

    crn_plt, crn_g=plot_c2rng()
    mrn_plt, mrn_g=plot_mrng(crn_g)

    rng_canvas = canvas.canvas()
    rng_canvas.insert(crn_plt)
    rng_canvas.insert(mrn_plt)
    rng_canvas.writetofile('mc-select-bias-range.eps')
    rng_canvas.writetofile('mc-select-bias-range.pdf')




do_all_plots()
