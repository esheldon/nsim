from __future__ import print_function

import biggles

thresh=[0,7,10,13,16,19]
s2n_means= [8.68, 15.47, 97.44]

xrng = [-2, 21]
log_xrng=[5, 140]

def add_mbounds(plt, xb, maxbias):
    mrng=biggles.FillBetween(xb, [maxbias]*2,
                        xb, [-maxbias]*2,
                        color='grey80')
    plt.add(mrng)

    z=[0,0]
    zline=biggles.Curve(xb, z)
    plt.add(zline)

def doplot(plt,xvals, yvals, yerr, yvals_corr, yerr_corr):

    nocorr_cur=biggles.Curve(xvals, yvals, color='red', type='solid')
    nocorr_pts=biggles.Points(xvals, yvals, color='red', type='filled circle')
    nocorr_err_pts=biggles.SymmetricErrorBarsY(xvals, yvals, yerr, color='red')
    nocorr_pts.label='uncorrected for selection'

    corr_cur=biggles.Curve(xvals, yvals_corr, color='blue', type='dashed')
    corr_pts=biggles.Points(xvals, yvals_corr, color='blue', type='filled circle')
    corr_err_pts=biggles.SymmetricErrorBarsY(xvals, yvals_corr, yerr_corr, color='blue')
    corr_pts.label='corrected for selection'

    key=biggles.PlotKey(
        0.9, 0.15,
        [nocorr_pts,corr_pts],
        halign='right',
    )

    plt.add(
        nocorr_pts, nocorr_err_pts, nocorr_cur,
        corr_pts, corr_err_pts,corr_cur,
        key,
    )

def plot_mrng():
    plt=biggles.FramedPlot()
    plt.xrange=log_xrng
    plt.yrange=[-8,8]

    plt.ylabel=r'$m [10^{-3}]$'
    plt.xlabel=r'$\langle S/N \rangle$'
    plt.xlog=True

    nocorr=[-0.37, 7.44, 0.47]
    err_nocorr=[0.55, 0.37, 0.23]

    corr=[-0.11, -0.41, -0.02]
    err_corr=[0.55, 0.37, 0.23]

    add_mbounds(plt, log_xrng, 1.0)
    doplot(
        plt,
        s2n_means,
        nocorr, err_nocorr,
        corr, err_corr,
    )

    plt.show()
    plt.write_eps('m-select-bias-range.eps')


def plot_mthresh():
    plt=biggles.FramedPlot()
    plt.xrange=xrng
    plt.yrange=[-4,4]
    plt.ylabel=r'$m [10^{-3}]$'
    plt.xlabel=r'$S/N$ threshold'

    mthresh_corr=[-0.91, -0.15, -0.10, -0.09, -0.04, -0.23]
    mthresh_err_corr=[0.19,0.19,0.20,0.20,0.21,0.22]

    mthresh_nocorr=[-0.91, 2.53, 3.48, 2.61, 1.85, 0.98]
    mthresh_err_nocorr=[0.19,0.19,0.20,0.20,0.21,0.22]

    add_mbounds(plt, xrng, 1.0)

    doplot(
        plt,
        thresh,
        mthresh_nocorr, mthresh_err_nocorr,
        mthresh_corr, mthresh_err_corr,
    )

    plt.show()
    plt.write_eps('m-select-bias-thresh.eps')

def plot_c2rng():
    plt=biggles.FramedPlot()
    plt.xrange=log_xrng
    plt.yrange=[-3,3]
    plt.ylabel=r'$c [10^{-4}]$'
    plt.xlabel=r'$\langle S/N \rangle$'
    plt.xlog=True

    nocorr     = [1.78, 1.16, 0.46]
    err_nocorr = [0.27, 0.18, 0.12]

    corr       = [0.27, 0.09, 0.20]
    err_corr   = [0.27, 0.18, 0.12]

    add_mbounds(plt, log_xrng, 1.0)

    doplot(
        plt,
        s2n_means,
        nocorr, err_nocorr,
        corr, err_corr,
    )

    plt.show()
    plt.write_eps('c-select-bias-range.eps')


def plot_c2thresh():
    plt=biggles.FramedPlot()
    plt.xrange=xrng
    plt.yrange=[-3,3]
    plt.ylabel=r'$c [10^{-4}]$'
    plt.xlabel=r'$S/N$ threshold'

    nocorr     = [0.17, 0.93, 0.81, 0.58, 0.49, 0.42]
    err_nocorr = [0.09, 0.10, 0.10, 0.10, 0.11, 0.11]

    corr       = [0.17, 0.17, 0.15, 0.13, 0.14, 0.12]
    err_corr   = [0.09, 0.10, 0.10, 0.10, 0.11, 0.11]

    add_mbounds(plt, xrng, 1.0)

    doplot(
        plt,
        thresh,
        nocorr, err_nocorr,
        corr, err_corr,
    )

    plt.show()
    plt.write_eps('c-select-bias-thresh.eps')


def plot_m():
    plot_mthresh()
    plot_mrng()

def plot_c2():
    plot_c2thresh()
    plot_c2rng()


plot_m()
plot_c2()
