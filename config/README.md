REAL METACAL
--------------

sim-bd01 and sim-bd01z
=======================
using galsim to do Gary's bulge+disk sim

- deep data
    - run-bd01zmcal-degrade01
        - shows significant e2, along direction of psf ellipticity
        - coellip3 psf modeling is probably shit
        - also the values are piling up at the boundaries in the two
            sided erf

    - run-bd01zmcal-degrade02
        - round the edges more in the two sided erfs
        - will address psf in next test
        - still additive error in e2

    - run-bd01zmcal-degrade03
        - round the edges more in the two sided erfs
        - em3 psf
        - still additive error in e2

    - run-bd01zmcal-degrade04
        - a bit of sub-pixel integration 4x4
        - still additive error in e2
        - no better: maybe it is all or nothing?

    - run-bd01zmcal-degrade05
        - calculating psf response

- regular runs
    - run-bd01mcal-t01
        - short run with frac err aimed to be ~0.0035
    - run-bd01mcal-t02
        meas: 0.0352182 +/- 0.000116463, 0.000519263 +/- 0.000115394
        fracdiff: 6.24e-03 +/- 3.33e-03

        additive is back

        sensitivity is much broader than bd02-04.  Either due to
        psf or due to centroid shift

        making cut s2n_r > 10 gives same answer!  must be a bug

    - run-bd01mcal-01
        - planned full run





- sim-bd02
- sim-bd02z
    - round gaussian psf, r50=1.5 same
    - no centroid shifts or bulge shift
- run-bd02zmcal-degrade01
    - rerunning after bug fix in sim
    - also bug fix in degradation: weight map was different
        for the regular and degrade run; noise was assumed to be 1/100 higher
        than it was.  could this make up the few parts in a thousand
        difference?
    - looks terrible, what happened?

- run-bd02mcal-t01
    * rerunning after bug fix

- so was it
    - the psf form
    - the psf ellipticity
    - the dev shift
    - the centroid shift

# test psf ellipticity
- sim-bd03
- sim-bd03z
    - elliptical gaussian psf, r50=1.5 same
    - no centroid shifts or bulge shift
- run-bd03zmcal-degrade01
- run-bd03mcal-t01
    - relatively short run
    - looks OK:  -0.0016 +/- 0.0022

# test dev shift
- sim-bd04
- sim-bd04z
    - elliptical gaussian psf, r50=1.5 same
    - add dev shift
    - no centroid shift
- run-bd04zmcal-degrade01
- run-bd04mcal-t01
    - relatively short run
    - this was it: bias of -9% as before
    - BUG FOUND in sim, not actually giving the galaxies shapes
        and shearing each individuall rather than the total
    -   0.03517 +/- 0.00012
       -0.00061 +/- 0.00012
    some detected additive bias
    fracdev on sheared component
        4.87e-03 +/- 3.28e-03
    1.5 sigma
- run-bd04mcal-t01b
    - more stats
        4.32e-03 +/- 2.38e-03
    - combined with t01
        4.51e-03 +/- 1.93e-03

- run-bd04zmcal-degrade02
- run-bd04mcal-t02
    - accidentally ran with prior g sigma 0.3 in both
        the degrade and normal run
    - I had modified it back but did not install it

        meas: 0.0352601 +/- 0.000116722, -0.000132118 +/- 0.000123185
        fracdiff: 7.43e-03 +/- 3.33e-03
    - additive seems to be gone now...

#  trying integration
- run-bd04zmcal-degint01
    - still additive errors same
    - turns out this was model bias, see below

- sim-bd05
    - high s/n
- run-bd05zmcal-degrade01
- run-bd05mcal-t01
    meas: 0.0349752 +/- 8.7483e-05, -7.44605e-05 +/- 8.38404e-05
    fracdiff: -7.08e-04 +/- 2.50e-03

    There is no detected additive or multiplicative error.
    So metacal is not correcting for noise properly 

    Maybe my degrading has a bug

        # for this high s/n, various measures agree between the
        # degraded and noisy run

        >>> t['mcal_g'].mean(axis=0) - deep['mcal_psf_sens'].mean(axis=0)
            array([  3.28291796e-02,  -6.99118029e-05])

        >>> t['mcal_g'].mean(axis=0) - t['mcal_psf_sens'].mean(axis=0)
            array([  3.28291706e-02,  -5.41409336e-05])

        >>> deep['mcal_g_sens'].mean(axis=0)
            array([[  9.38643015e-01,   8.87911118e-04],
               [  4.66034495e-06,   9.41100525e-01]])

        >>> t['mcal_g_sens'].mean(axis=0)
            array([[  9.37297464e-01,   5.92893438e-04],
               [ -3.01999789e-04,   9.42193738e-01]])

        >>> t['g'].mean(axis=0)
            array([ 0.03275939,  0.00289699])

        >>> t['mcal_g'].mean(axis=0) - deep['mcal_psf_sens'].mean(axis=0)
            array([  3.28291796e-02,  -6.99118029e-05])

- lets simplify
- sim-ggnr02
- sim-ggnr02z
    - gauss model and psf
- run-ggnr02zmcal-degrade01
    - at bnl
    - Rpsf is well detected but small
            -8.06407392e-05 +/- 4.20040477e-06

- run-ggnr02mcal-t01
    - at slac
    - additive is huge (and Rpsf is measured to be even larger)
        meas: 0.0351855 +/- 0.000108188, 0.00190488 +/- 0.000110891
        fracdiff: 5.30e-03 +/- 3.09e-03

    so something is going wrong in the degradation.  note in some runs
    the degradation (even for skynoise/1000 start) was at least working
    reasonably well in reproducing the Rpsf

    - ideas
        - need larger start noise, to include some of the issues?  In
            my other sims I used noise/10
            - looks the same in a quick test
        - maybe difference running at bnl vs slac
            - the degrade ran at bnl and the regular at slac.  So try
                running the degrade again at slac to see if suddenly the
                additive term appears.

- run-ggnr02zmcal-degrade02
    - at slac to see if additive appears
    - quick test shows (on 1000 galaxies) that it is not there, or
        at least small 8.69472673e-05 +/- 3.74688387e-05
    - trying a deepcopy on obs in case I missed something
        - no

- ideas
    * look at high s/n
    * maybe the prior on g is too restrictive.  There is a population of high
        ellip objects when the bulge is offset
    * the pixel integration might still help at the level of 4 parts in a thousand.
        Should check on a high s/n run
    * centroid might move somewhat differently for the degraded.  Could either
        loosen prior or tighten it to see what happens
    * sim might want to keep center of light near center of image

exploring the psf leakage

- sim-ggnr01
    - psf 0.0,0.05
- run-ggnr01max01
    - high s/n fitting correct model
    - see if additive bias
    - looks OK

- sim-egnr06
    - psf 0.0,0.05
- run-egnr06max01
    - high s/n
    - fit gauss to exp
    - see if additive bias
        >>> t['g'].mean(axis=0)
            array([ 0.07516843,  0.00255811])
        >>> t['g'].std(axis=0)/sqrt(t.size)
            array([ 0.00029705,  0.00029721])
      !! that was it:  model bias can lead to psf leakage!!

sim-eg17
===========

- run-eg17mcal07
    - s/n=100
    - straight
         0.0035 +/- 0.0066

- run-eg17mcal08
    - s/n=100
    - from dg06mcal01 rerun, wonder if something has gone wrong
    - also more stats
    - straight average
        0.00078 +/- 0.0031

- run-eg17mcal09
    - fitting gauss
    - s/n=100
    - straight
        -0.0013 +/- 0.0031
    - reran
        - g_mean as field
        -straight
            0.0038 +/- 0.0033
        -weighted
            -0.0011 +/- 0.0032



- run-eg17mcal10
    - gauss fit
    - s/n=25
    - straight
        0.054 +/- 0.0031
    - re-running with full rendering
    - straight
        0.049 +/- 0.0033
    - weights
        0.031 +/- 0.0031
    - rerunning with new code, but I don't expect change
        -straight
            0.051 +/- 0.0033
        - weights, g_mean
            0.04 +/- 0.0031
    - big jump from s/n=50

- run-eg17mcalt02
    - exp fit
    - full rendering
    - s/n=25
    - straight
        0.048 +/- 0.0032
    - weighted
         0.037 +/- 0.0031

    - very similar to when fitting a gaussian. I think this indicates the
        problem is just related to the noise not being identical for the
        different images

- run-eg17mcalt03
    - exp fit
    - whitening

- run-eg17mcalt05
    - exp fit
    - s/n=25
    - just resimulating the model.  should only correct noise bias,
        so need to run with true model
    - undercorrecting from quick test
    - forgot will be *very* biased location
        -0.039

- run-eg17mcalt06
    - added lots of outputs to calculate e.g. gmean,
        gnoshear, sens from model
    - see that gmean*0.5*(1/sens + 1/sens_model) is nearly unbiased
- run-eg17mcalt07
    - same os 06 but fitting gauss
    - doesn't work

- run-eg17mcalt08
    - using moms
    - crap at high s/n

- run-eg17mcalt09
    - using em
    - crap at high s/n

- run-eg17mcalt10
    - using admom
    - for deconv. psf version, works at high s/n to 0.25%
        - would only work for gaussian psf
    - just fitting a gaussian convolved with gaussian psf works
        better

- run-eg17mcalt11
    - gauss fit
    - s/n=50
        - weights -0.00065 +/- 0.0038
        - straight 0.0099 +/- 0.004

- run-eg17mcalt12
    - gauss fit
    - s/n=10000
    - to be used for nearest

- run-eg17mcalt13
    - gauss fit
    - nearest from t12 above
    - s/n=100
        - works
    - s/n=50
        - straight
            0.0064
        - weights
            -0.0034
    - rerunning with BA
        - s/n=100
            0.0017 +/- 0.0026 (6% correction)
        - s/n=50
            -0.003

- run-eg17mcalt14
    - gauss fit
    - nearest from t12 above
    - s/n=23 and higher
        - decent with weights
        - one big outlier... similar to what I saw before. Maybe
            remove huge outliers
    - rerunning with BA

-run-eg17zmcal01
    - zero shear training run

- run-eg17mcalt15
    - gauss fit
    - nearest from zero shear eg17zmcal01
    - many s/n
    - wierd s/n=23 for BA, does't fit with
        others; more outlier stuff?

** gaussmomn
- run-eg17zgmom01
    - gauss mom fitting, high s/n for training
- run-eg17gmomt01
    - using run-eg17zgmom01 for training.  Realized this
        full p(g) will necessarily be "no ring". switching
        to sim-egnr05

- run-egnr05gmomt01
    - using run-eg17zgmom01 for training

sim-gg11 and sim-gg1z
=====================

For testing moment fitting and prior template sampling

- run-gg11zgmom01
    - high s/n for training
    - aiming for 100,000 templates

- run-gg11gmomt01
    - expanding true shear, h=1.0e-3
    - mean shear comes out about right (zero before adding truth
        back in) but the covariance matrix is crazy.  Maybe
        it was zero just due to errors

- run-gg11gmomt02
    - don't expand true shear, this should not come out zero in
        raw shear
    - seeing shear in component 2 as well; a clue!
- run-gg11gmomt03
    - h=1.0e-6
    - nope, same

- run-gg11gmomt04
    - s/n=10 even though gaussian assumption breaks down
    - should give stable results
    - errors look reasonable, but seeing signal in e2
        - was bug in sim

* run-gg11gmomt05
    - s/n=12, expand true, still draw truth for templates
    - no ring
    - still terrible, maybe it is draw truth
        - re-running without seed bug
        - still terrible

- run-gg11gmomt06
    - s/n=12, expand [0,0], draw truth for templates, 
        except centers which are randomized uniformly
    - no ring
    - also improves stability of isampling a bit, but still
        getting low neff
    - looks OK, still need
        - better sampling
        - proper likelihood functions

- run-gg11gmomt07
    - same as run-gg11gmomt06 but using the training
        set run, nrand_cen=10
    - OK but seeing some significant shear 2!
        - was bug in sim

- h 1.0e-3 maybe use 1.0e-6?


sim-dg06
===========

- run-dg06mcal01
    - s/n=100
    - fit dev
    - straight
        -0.0046 +/- 0.0051
    - weights
        -0.0037 +/- 0.0051
    - did new run with more stats
    - straight
        0.0054 +/- 0.0025
    - that's a big swing... something wrong?

- run-dg06mcal02
    - s/n=100
    - fit gauss
    - ~1.4% bias
    - putting sensitivity range [0.0, 2.0] results
        in 0.0047 +/- 0.0058
    - adding weights with SN=0.22 0.0092 +/- 0.0059
    - both 0.003 +/- 0.0057

- run-dg06mcal03
    - s/n=100
    - fit exp
    - straight
        0.0067 +/- 0.0027
    - weights SN=0.22
        0.0054 +/- 0.0027
    - sens [0.0, 2.0] and weights
        0.0022 +/- 0.0027

- run-dg06mcal04
    - s/n=50
    - fit exp
    - straight
        0.027 +/- 0.0022
    - weights SN=0.22
        0.0054 +/- 0.0027
    - sens [-0.5, 2.0]
        0.019 +/- 0.0022
    - sens [-0.25, 2.0] and weights
        0.996 kept
        0.0044 +/- 0.0021
    - sens [-0.15, 2.0] and weights
        0.993 kept
        -0.0061 +/- 0.0022
        
    - sens [0.0, 2.0] and weights
        0.95 kept
        -0.032 +/- 0.0021
    - sens [0.0, 2.0] no weights
        same

- run-dg06mcal06
    - exp
    - very high s/n 10000
    - straight average
        9.2e-05 +/- 0.00014
    - sensitivity peaks at abot 0.8!

- run-dg06mcalt01
    - test run with exp fit
    - s/n=100
    - using full rendering of models
    - straight
        0.0047 +/- 0.0027
    - weighted
        0.0033 +/- 0.0027
- run-dg06mcalt02
    - fit exp
    - s/n=50
    - using full rendering of models
    - straight
        0.023 +/- 0.0026
    - with weights
        0.021 +/- 0.0027

- run-dg06mcalt03
    - fit exp
    - s/n=50
    - back to trying unsheared, just to see difference with 02
        0.029 +/- 0.0027

- run-dg06mcalt04
    - fit exp
    - isample
    - s/n=50
    - straight
        0.022 +/- 0.0028

- run-dg06mcalt04
    - fit exp
    - isample
    - s/n=50
    - using original fit shape
    - straight
        0.013 +/- 0.0027
    - weights
        0.011 +/- 0.0027

# using prior info
- sim-dg06z
    - zero shear, calculating metacal pars
    - will use for priors
- run-dg06zmcal01
    - fitting exp
    - will use for priors

- run-dg06mcalt06
    - fitting exp
    - bring in prior from run-dg06zmcal01
    - s/n=100 fine
    - s/n=50 fine
    - running again to make sure nothing changed
        weights
        fracdiff(lensfit): -0.005 +/- 0.0026
        fracdiff(pqr):     -0.0023 +/- 0.0022
        no weights
        fracdiff(lensfit): -0.0028 +/- 0.0027
        fracdiff(pqr):     -0.0024 +/- 0.0022



- run-dg06mcal07
    - fitting exp
    - bring in prior from run-dg06zmcal01
     [ 23, 35, 53, 81, 123, 187, 285, 433, 658, 1000 ]
    - 23 blew up, why?

- run-dg06mcal07-b
    - ran with larger allowed errors, still blows up
        at s/n=23

- run-dg06mcal08
    - fitting exp
    - bring in run-dg06zmcal01
    - using sensitivity directly in lensfit sum
    - 100 looks fine
    - implemented prior from the run.
    - rerunning at s/n=100 to make sure all went well
        - ok 
               weights: -0.0016 +/- 0.0055
            no weights: -0.00068 +/- 0.0056
    - now running at s/n=50
        - lensfit: -0.0079 +/- 0.0053
        - BA using mean response: -0.0042 +/- 0.0053

- run-dg06mcal09
    - same as run-dg06mcal08 but at s/n=23
    - note pqr does not have any corr intrinsic to it,
        so we do --corr-model to apply a correction
    - pqr with model bias post-fix
        0.0033 +/- 0.0053
        - maybe using correct priors was what was needed
    - lensfit with model bias included in sens. calc.
        -0.023 +/- 0.0051
    - may want to do a run where we keep both lensfit
        style sensitivities, and model bias

- run-dg06mcal10
    - same as 09 but keeping both kinds of lensfit sens
    - pqr post model corr
        0.011 +/- 0.0042
        wtf?  sometimes get those fluctuations
    - lensfit with g_sens_r
        -0.0147
    - lensfit with post model corr
        no weights: 0.017 +/- 0.0053
        weights: 0.017 +/- 0.0052

- run-dg06mcal11
    - doing isampler iterattion andupping the number of isamples, saving
        neff [500,4000]
    - s/n=23 same as 09 and 10
        pqr model corr post
            0.0044 +/- 0.0053
        lensfit model corr pose
            0.01 +/- 0.0053
        lensfit g_sens_r
            0.0098

- run-dg06mcal12
    - same as run-dg06mcal11 but s/n=10
    - lots of jobs died it seems
    - corr model post (note no jack on pqr, so errors too small)
        fracdiff(lensfit): -0.0024 +/- 0.0064
        fracdiff(pqr):     -0.0087 +/- 0.0038

        adding back in the jobs that died
        fracdiff(lensfit): 0.0017 +/- 0.0053
        fracdiff(pqr):     -0.0045 +/- 0.0031

        lensfit with weights: 0.0082 +/- 0.0053

        The errors for pqr are also probably about 0.0053 when jackknifed. So
        the errors are probably a bit underestimated, maybe because there is
        significant spread in this mean model correction; it is hard to say.

        Maybe errors being off is due to using the convolved image?

        But it looks decent.  So the detail I needed to get right was the
        prior being from the same model we are using at low s/n. Using
        more samples may also have been important.

        Need to get pqr jackknifing working from the chunks

* proper weighting for pqr etc.
* can do "importance sampling" over the template set,
    instead of some other likelihood/posterior sampling

- run-dg06mcalt07
    - fit gauss
    - s/n=50

- run-dg06mcalt08
    - metacal-psamp
    - fit exp
    - s/n=50
    - new prior sampling technique, sampling from actual
        parameters measured from a high s/n run

- run-dg06mcalt10
    - metacal-isample-nearest
    - fit exp
    - s/n=50
    - priors from the high s/n fit
    - may want to remove the prior we apply, but then calculate
        sensitivity from the true prior (or best guess in real
        data)?
    - 0.004 in both
    - running again (on slac)
        fracdiff(lensfit): 0.0074 +/- 0.0026
        fracdiff(pqr):     0.0059 +/- 0.0022

- run-dg06mcalt11
    - metacal-isample-nearest
    - fit gauss
    - s/n=50
    - priors from the high s/n fit
    - may want to remove the prior we apply, but then calculate
        sensitivity from the true prior (or best guess in real
        data)?

        fracdiff(lensfit): 0.0085 +/- 0.0028
        fracdiff(pqr):     0.0081 +/- 0.0024

- run-dg06mcalt12
    - metacal-isample-nearest
    - including centroid in nearest match
    - true priors
    - fit exp
    - s/n=100
        fracdiff(lensfit): 0.0021 +/- 0.0027
        fracdiff(pqr):     0.0019 +/- 0.0022
    - s/n=50



- sim-dg07z
    - zero shear, no ring sim used as metacal training for sim-dg07
    - run-dg07zmcal01
        - using still max like, plus adding Richardson extrapolation
    - do a non-max like run?  Might matter at the part in 1000 level

- sim-dg07
    - shear = 0.03 for testing metacal.
    - We expect possible bias at
        the 0.0015 level from metacal failing at higher shear
    - expect B&A to have error also, at 0.002 level
    - both are postive, so we expect we need to expand B&A around
        true shear to avoid that bias

    - at higher shear the nearest neighbor matching won't work as well
        - argument for doing this with full sheared version of the
            prior rather than just using derivatives.

- run-dg07mcaltest01
    - at s/n = 100 expanding about true shear.
        0.0051 +/- 0.0018
    - using mean sensitivity from templates 0.004, so maybe
        the space is too sparse?
    - s/n=50
        - pqr
            - like averaged 0.0030 +/- 0.0018
            - overall mean correction 0.0028
    - another run using new zero shear training set dg07z
        0.0035 +/- 0.0018
      using total mean from training set
        0.0026

- run-dg07mcaltest02
    - pqr
        0.0074 +/- 0.0014

* match pars in regular fit instead of convolved fit?
* 
* use mean relation for response
    - not clear it will help much

- new idea: just degrade high s/n images using same noise image
- templates need to be done that match the data set at hand.
    For my sims I'll just redo for each s/n level

- s/n=50
    - run-dg07zmcal-degrade50
        - run with s/n=500 degraded to 50
    - run-dg06mcal50-01
        - run at s/n=50 for matching above using g_noshear
            -0.00127 +/- 0.00531
        - using g_mean
            -0.0035 +/- 0.0053
        ( using just g I get 0.00044 but that's probably a fluke )
    - run-dg06mcal50-02
        - more precision
            0.00118 +/- 0.00106

- s/n=20
    - run-dg07zmcal-degrade20
        - run with s/n=200 degraded to 20
        - hmm... do we need more samples to get the mean correction
            precise enough?
                0.65859 +/- 0.00052

    - run-dg06mcal20-01
        - run at s/n=20 for matching above using g_noshear
            - using g_noshear
                 0.0045 +/- 0.0053
            - using g_mean
                -0.0013 +/- 0.0053

    - run-dg06mcal20-02
        - more precision
            - g_noshear
                0.0032 +/- 0.0021
            - g_mean
                -0.0028 +/- 0.0021

- s/n=10
    - run-dg07zmcal-degrade10
        - run with s/n=200 degraded to 10

    - run-dg06mcal10-01
        - quick run, expect error in fracdiff of 0.005
        - using g_noshear
        - using g_mean
            0.0030 +/- 0.0053

    - run-dg06mcal10-02
        - more precision
        - g_mean
            -0.00027 +/- 0.00214
        - g_noshear
            0.0090
        - I wonder if this points to something cancelling.  Also
            maybe ring test biting me again?

higher shear
-------------

sim-dg01,sim-dg08
==================
same as sim-dg06 but shear=0.08,0.05 respectively
still run-dg07zmcal-degrade50 for training set

- run-dg08mcal50-01
    - shear 0.05
    - s/n=50 lots of stats

        - using g_mean
           -0.000979 +/- 0.000857


- run-dg01mcal50-01
    - shear 0.08
    - s/n=50 lots of stats

        - using g_mean
            0.0012 +/- 0.00054

- run-dg01mcal50-02
    - for more stats
        0.00021 +/- 0.00055
- combined dg01mcal50-01 and 02
    0.000721 +/- 0.000383
    7.21e-4 +/- 3.83e-4

sim-dg09
=========
shear [0.035,0.035]
- run-dg09mcal50-01
    - now doing sensitivity for both components
- run-dg09mcal50-02
    0.00063 +/- 0.00061  0.00094 +/- 0.00061

- combined -01 and -02
    0.00084 +/- 0.00054  0.00053 +/- 0.00054
    combined using more stats training set "run-dg07zmcal-degrade50-lots"
    5.4e-04 +/- 5.4e-04  3.0e-04 +/- 5.4e-04

- run-dg07zmcal-degrade10-lots
    - on comet
- run-dg09mcal10-01
    - not yet run

- is it the T/F or g priors?

    - "true" priors applied for deep run, although they may not have
        mattered much
        - the actual recovered values very different from priors
    - for runs that worked same priors were applied.

    - the sensitivities are coming out very low compared to the
        run with true g prior
        - model sens less different
        - should I apply true prior?  Then the issue is I have
            no examples with which to calculate the model sensitivity in
            regions where the true prior differs
    - this is the problem with using biased measurements

    - I didn't match prior galaxies in centroid












s/n=50 results with no cuts
---------------------------

sim-egtight01: broad distribution in cen, not needed
sim-egtight04: tight distribution on cen, probably the right thing to have always

run              |   priors          | g during?          | max_fitter | nsamp | df  | fracerr              | comments
-----------------+-------------------+--------------------+------------+-------+-----+----------------------+---------
egtight01isamp04 | flat priors       |  g prior during    | lm         | 2000  | 2.1 |  0.016 +/- 0.0014    |
egtight01isamp06 | flat priors       |  g prior during    | lm         | 4000  | 2.1 |  0.018 +/- 0.0014    |
egtight01isamp08 | flat priors       |  no g prior during | lm         | 2000  | 2.1 |  0.017 +/- 0.00071   | note smaller error
egtight01isamp09 | flat priors       |  no g prior during | lm         | 4000  | 2.1 |  0.017 +/- 0.00066   | note smaller error
egtight01isamp11 | flat priors       |  g prior during    | lm         | 24000 | 1.0 |  0.017 +/- 0.0014    |
egtight01isamp12 | true TF priors    |  g prior during    | nm         | 2000  | 2.1 |  0.0072 +/- 0.00071  |
egtight01isamp13 | true TF priors    |  g prior during    | nm         | 2000  | 1.0 |  0.0078 +/- 0.0007   |
egtight01isamp14 | true TF priors    |  g prior during    | nm         | 4000  | 1.0 |  0.0067 +/- 0.00068  |
egtight01isamp15 | true TF priors    |  g prior during    | lm         | 2000  | 1.0 |  0.0066 +/- 0.00068  |
egtight01isamp16 | fixed flat priors |  no g prior during | nm         | 2000  | 2.1 |  0.017 +/- 0.0007    | right priors for log e
egtight01isamp17 | fixed flat priors |  no g prior during | lm         | 2000  | 1.0 |  0.016 +/- 0.00069   |
egtight01isamp18 | fixed flat priors |  no g prior during | nm         | 2000  | 1.0 |  0.015 +/- 0.00069   | uninformative prior on cen this time

egtight04isamp01 | flat/weak priors  | no g prior during   | nm        | 2000 |  2.1 | 0.016 +/- 0.00071    | no better!
egtight04isamp02 | true TF priors    | no g prior during   | nm        | 2000 |  2.1 | 0.0059 +/- 0.0007    | note noise error is much larger than the cen distribution
egtight04isamp03 | true TF priors    | g prior during      | nm        | 2000 |  2.1 | 0.0072 +/- 0.00068   |
egtight04isamp04 | flat/weak priors  | g prior during      | nm        | 2000 |  2.1 |                      |

egtight04isamp06 | flat/weak priors  | g prior during      | nm        | 2000 |  2.1 |                      | limit errors more

egtight04isamp07 | true TF priors    | no g prior during   | nm        | 2000 |  2.1 | 0.0058 +/- 0.0007    | draw truth: as expected, same result
egtight04isamp08 | true TF priors    | no g prior during   | nm        | 2000 |  2.1 |                      | draw truth and set isample weights to 1. Should be unbiased

egtight04mc01    | true TF priors    | g prior during      | nm        | -    |  -   | 0.0058 +/- 0.0071    | mcmc fitter, consistent with the isampler

broader distributions
eg04isamp01      | true TF priors    | g prior during      | nm        | 2000 |  2.1 | 0.0086 +/- 0.00077   | with weights -0.0077 +/- 0.00072

s/n=25
egtight04isamp05 | flat/weak priors  | g prior during      | nm        | 2000 |  2.1 | 0.041 +/- 0.00073  | with weights 0.031 +/- 0.00073
egtight04isamp09 | true TF priors    | g prior during      | nm        | 2000 |  2.1 | 0.013 +/- 0.00074  | with weights -0.00052 +/- 0.00073

broad
eg04isamp02      | true TF priors    | g prior during      | nm        | 2000 |  2.1 | 0.013 +/- 0.00079  | with weights 0.0053 +/- 0.00076

eg04isamp03  | 12 s/n bins, ifactor=1.414, df=2.1, lm, all true priors       | 1.5% bias at s/n=10
eg01isamp01  | 12 s/n bins, ifactor=1.414, df=2.1, lm, all true priors       | looks good
eg01isamp03  | 12 s/n bins, ifactor=1.000, df=2.1, lm, flat on flux          | looks like crap
eg01isamp04  | 12 s/n bins, ifactor=1.414, df=2.1, lm, flat on flux          | looks a little better (note weights do not help)
eg01isamp05  | 12 s/n bins, ifactor=1.414, df=2.1, lm fix cov, true priors   | looks good
eg01isamp06  | 12 s/n bins, ifactor=1.000, df=2.1, lm fix cov, true priors   | increases efficiency, bias about same
eg01isamp07  | 12 s/n bins, ifactor=1.000, df=1.0, lm fix cov, true priors   | lowers fracuse (argh, meant for flat test) maybe a bit worse?

s/n=50
eg01isamp08  | s/n=50, ifactor=1.000, df=2.1, lm fix cov, true priors, *iter* | no eff. improvement for sim
eg04isamp04  | s/n=50, ifactor=1.000, df=2.1, lm fix cov, true priors, *iter* | no eff. improvement for sim


low nsample
eg01isamp09  | nsample 500                 | looks fine!
eg04isamp05  | nsample 500, ifactor=1      | looks different from eg04isamp03, is it samples or ifactor?
             |                             | also used nm, larger err max for isamp03
eg04isamp06  | nsample 500, ifactor=1.414  | kind of in between

eg01isamp11  | 12 s/n bins, ifactor=1.000, df=2.1, lm fix cov, true priors   | recovering from bug

eg04isampt   | s/n=50, ifactor=1.414, df=2.1, iter 500,2000, true priors, prior *not during* | looks about as expected

new sampling on grid of p
-------------------------
eg01isampp01 | nsample 500 | 
egnr01isampp01 | s/n=1000 no ring, nsample 500 | agrees with lensfit
egnr01isampp02 | s/n=50 no ring, nsample 500 | 
egnr01isampp03 | 12 s/n values, 500 |

egnr02isampp01 | 12 s/n values, 2000 |
egnr02isampp02 | 2-d grid, s/n=1000 |

#
# no ring
#

egnr02isamp01  | s/n=50, ifactor=1.000, df=2.1, lm fix cov, true priors | testing s/n bias
egnr02isamp02  | s/n=50, ifactor=1.000, df=2.1, lm fix cov, true priors | added saving round

exp + dev
---------
sim-egdg01

egdg01isamp01 | 12 s/n bins, "true priors" but only fitting exp | significant bias
egdg01isamp02 | 12 s/n bins, "true priors" fitting composite model | save psf s/n


GPriorMErf sim


# todo

* run egcosmos and BA prior
* Figure out more realistic T and F distributioins for sim.  Might be
  able to get away with lognormal with very large variance and
  where the *peak* is at say T=4 and F=whatever.  but the sim
  would have to be updated for that, since it assumes the mean
  for s/n
* might want to do an intial max fit in linear space: if T is <= 0
    then don't continue.  Or at least save that for later.
- test with linear F: sucks
- test with linear T: sucks
* get isamp working for B&A: weights
* run on great-des
* run on testbed

# trying log T only
    - changes in _gmix.c
    - in set_prior
    - in the guess for maxlike for isample
    * note I think the guessers were not doing log for feeding
        into maxlike for emcee....

test s/n=1000 for 
    - BA pdf and BA prior
        + eg04rtest07
    - cosmos pdf and BA prior
        * egcosmos01rtest02
    - cosmos pdf and cosmos prior (seems to have instability)
        * egcosmos01rtest03

fwhm=1.2
    - sim-eg04
- new Observation and log style
    - sim-eg12
    - run-eg12rtest1
        - s/n=100

- great3 priors
    - look good as long as I keep the g prior separate.  Calling
    this "hybrid"

    - sim-eg11
        - exp fits from the real_galaxy deep as prior, things look good.  
    - todo is to go back to bdf and do the g_prior separately and see if I can
      do decently well.  Will want to fit bdf on the real_galaxy deep fields.


* general comments so far
    - get basically same results for draw priors, truth, maxlike

    - I see a "wave" pattern as a function of s/n. I think this was bad seeds.
      Turns out mtrand uses /dev/urandom.  Had to replace emcee._random with my
      own.

        - gg04rtest
            - using guess around the truth to see if the wave disappears.  If
              so we can try a max likelihood guess?
                    - s/n=81 looks the same
            - re-used gg04rtest for a=2. a little bit better, suggestive.
              2-sigma from truth might warrant a new run to bring up s/n

        - gg04rtest3 same configuratioin as gg04rtest

            - a=2? a little bit better, suggestive. 2-sigma from truth
               - might warrant a new run to bring up s/n
            - more steps or walkers?
            - Try MH?  I have no indication this would work better
                - gg04rtest2
            - try isample? I have no indication this would work better

- coellip4
    - test on exp profile
        - eg01rctest1 guessing with broad range around truth; 0.05 scatter in
          the ellipticity.  s/n=50.  fracdiff = -2.18e-04 +/- 7.56e-04

        - eg01rctest2 guessing based on an exponential LM fit to the profile,
          which itself has a biased guess.  This is I think a much more stringent
          test, since I have seen the lm fit can be quite far off.
            - terrible, percent level

    - test on dev profile
        - do LM guess using exp as well...

current set
    sigrat2.0 gg01rcomb eg01r01 dg01rcomb
    sigrat1.4 gg03r01   eg03r01 dg04rcomb
    sigrat1.0 gg04rcomb eg04r01 dg05r01


J have jackknife errors calculated
- sigrat 2
    - sim-gg01
        - gg01r01  usual 80,400,200 prior after
            - looks good, although error bars seem too small
        - gg01r02 prior during
            - looks about the same
        - gg01rtest a=2 and lots of stats; looks pretty good at fracdiff -2.8e-4
          +/- 1.25e-4.  Could live with that!

        - J gg01r03 same as the test run, everything same as gg01r01 but with a=2
          also smaller error by factor of 2/5. Looks quite good; do see the "wave"
          at low amplitude.
        - J gg01r04 guess from truth, a=2, running
        - J gg01rcomb average of r03 r04

    - sim-eg01 sigrat 2
        - J eg01r01
        - looks quite good.  Problems at high s/n maybe from slow burnin?

    - sim-dg01 sigrat 2
        - J run-dg01r01
        - J run-dg01r02 guess true, a=2, running
            - errors are larger, but looks rather similar

        - J dg01rcomb combined 01 and 02

        * dg01r03 new seeding

- sigrat 1.4
    - sim-gg03
        - J run-gg03r01
    - sim-eg03
        * J run-eg03r01
    - sim-dg04
        - J run-dg04r01
        * J run-dg04r02 redoing with new seeding
        * J run-dg04r03 with new seeding
            - averaging...
        - J run-dg04rcomb currently holds old average
        - J run-dg04rcomb2 currently holds new average from all 3

           

- sigrat 1
    - sim-gg04
        - J gg04r01
            - strong "wave" pattern.
        - J gg04r02
            - try with new settings draw truth and a=2, still high error 5.0e-5
            - looks pretty consistent with gg04r01
        - J gg04rcomb averaged those
    - sim-eg04
        - J run-eg04r01
            - similar wave pattern
    - sim-dg05
        - J run-dg05r01, running using maxlike
            - pretty bad at low s/n

        - run-dg05rtest1 s/n 15 only  ..!! looks fine!
        - run-dg05rtest2 same         .. looks bad
        - run-dg05rtest3 maybe it is the deriviatives.  Try prior during to
          smooth it out?  Also might try to resurrect the mathematica stuff to
          get analytic
        - run-dg05rtest4 now have analytic deriviatives at non-zero shear! 4
          times slower than the numerical ones, but this can be optimized

        - pretty clear the errors are underestimated.  Will re-set the
          predicted error bars based on jackknife

        - J run-dg05r02. running. Predicted errors now from jackknifed
          the run-dg05r01 
            - looks similar
        
        - J run-dg05r03 new seed from /dev/random, some minor changes to priors
          functions.  Going a bit larger error, 7.5e-5
            - ARGHHH turns out emcee will start it's own random number generator
            so this did nothing
        - J run-dg05r04
            - same as 03 but now actually seeding from /dev/random
            - certainly looks random but errors are large

FWHM RAT
---------

- FWHM 1.2
    - turns out for dev srat 2 is same as fwhm rat 1.2, very far off for others

    - sim-eg05 0.08
        - eg05rtest1
        - srat 0.85, T rat of 0.715.

        - hypothesis: higher ellipticity is having problems due to an aperture
          effect.  Would be worse for higher shear. Why ok at higher s/n?

        - hypothesis: psf sampling is becoming a problem?  Why OK at higher s/n?

        - hypothesis: there are too many small ones to be realistic. 
            - trying again with smaller width for both flux and size
            - shear 0.08
    
    - sim-eg06 shear 0.15
        - did mean over the ellip distribution, 
        - might be nsigma render
        - use 10 instead of 5
        - nope same
        
        - trying larger psf
            - looks a bit different, higher instead of lower.  Maybe better?

    - sim-eg07 0.15
        - now psf T=30, sigma 3.87 to see if even more improvement
    - sim-eg08 0.08
        - recenter psf probably not right
        - no difference

        - try lower shear?  How to justify?
            - running

    - current picture
        - in convolved FWHM ratio space the exp seems to perform worse than dev
          at 1.2
            - this fwhm ratio corresponds to smaller exp in terms of moments
              than dev
            - but I haven't checked s/n=5 for dev
        - even with tight priors on size and flux biased at s/n=5
        - is prior too tight?  Would s/n flux or size be lower without it and
          indicate a problem?  Yes.  Check cosmos for width of prior?
        - still seeing waves in some runs but not others

flux limit based on cosmos
--------------------------

- set the flux mode to have s/n=5
- then just make size cuts
    - do in fwhm ratio space.  Means different cuts for exp and dev.
    - can also change s/n scale but need to think how best to do that, because
      the distribution is joint

older stuff
-----------

- MH
    - 10000 burnin, 10000 step
        - was more biased than emcee
        - gg01rtest6 did prior after and looks better, but still more biased
        - gg01rtest7 prior after and now 40000 after burnin to see if looks
        better
            - yeah, maybe a bit better.  So maybe it is really just need of
            more points?

- very high emcee steps            
    - reused run-gg01rtest to try a much larger number of points
    - nstep 2000x80!!  Same burnin 400x80
        - same
    -  maybe the gain was all in prior after.  Try fewer walkers=20
        run-gg01rtest2  quick look at plots doesn't look good

- tighter distribution in flux 
    - use 10% scatter instead of 30%.
        - sim-gg10
        - looks about the same

- try nearly-fixed other paramters besides shape
    - sim-gg02.  Widths 1.0e-5 for cen and frac for T and counts
    - gg02r01 looks crappy!  Is it because we only used 200 step? Doubt it
      since gg01r01 was fine. Simply number of trials? 1,700,000 should be enough..

      This is an important clue.  Could it be the sampler itself? Or how priors
      are calculated?

    - remake sim to use more reasonable distributions. Say cen_sigma 0.01 and
      width frac 0.01 for T and counts
      - run-gg02rtest1 with s/n=10


- fixcen
    - sim-gg09
    - run-gg01rtest1
        - s/n 10  80,400,200
        - looks the same pretty much

- isample
    - gg01rtest4
        - 200 min had bias high. Most were > 1000 actually
        - try 1000 min, same, maybe no surprise
        - try 5000 min same
        - try 200 but no prior during. Looks better. I wonder if I changed
          something else at the same time that made the difference.
          
          Re-running again because I *did* make a few changes after starting
          that run that were supposed to be irrelevant.  Looks like shit now.

          argh... now very low. running again rtest3 same pars to see if we
          are consistent. I don't know, could be OK

        - running s/n = 15 to see

- bdf
    - bdf01r01 terrible
        - was not keeping low arate, and many had low arate.
    - bdf01r02
        - better guessing and 2 pass.
        - looks much better but not great
    - bdf01r03
        - trying draw_truth to see if improved.  a=2, which might be an issue?
          Also prior during.

          ? try making prior *very* thin for dev side and then fitting
          both dev and b+d and taking best?

    - bdf01rtest re-starting at guess from best likelihood helps.
        - fewer with low arate
        - Unfortunately I also was using prior during so not absolutely sure
          which was most helpful.
    - bdf01rtest3 -
        - trying true guess with mca_a=2 first 20 pairs look good.
        - try with the "retry" version and not true guess.
            - some failures
            - trying a=3 first pass a=2 second, not enough
            - tried putting a chunk of bfrac near 0 and near 1 and it didn't
              fail once in 20 pairs.
            - trying a full run, but note we should test a few things.

        - try drawn priors and not forcing restart and a=3 and not removing low
          arate
            - looks terrible, 1.4% at s/n=35
        - trying nwalkers=320, a=2
            also rtest2 at s/n=23

    - also not was not keeping ones that failed arate, might want to change
      that?

- current runs
    - testing why small objects are biased
        - run-gg08rtest shear 0.01 and no expand shear true, test if was
          derivatives.  Might not be that useful until we do many such runs,
          since fractional error will be 1.e-3

- sim-dg01
    - run-dg01r33 shear 0.01 with nwalkers=20 to verify we still see high bias
        - I do see the bias!  So it was not just the better guess

- sim-dg02
    - run-dg02r01 nwalkers=40 and few s2n vals (should have used 1000 at the
      high end!)
    - run-dg02r02 same as r01 for more statistics

    - + run-dg02r03 using log spaced s/n value set

    - run-dg02rtest this can be re-used for other purposes
        - already saw going to 80 walkers does't help that much
        - try nwalkers=40 and mca_a=2?
            - looks like crap!  huh...

- sim-dg03  shear 0.08, see if OK at high shear
    - run-dg03r01 nwalkers=40, s/n starting at 10

    - + run-dg03r02 nwalkers=80
        - looks better!  So is it the nwalkers, burnin or the nstep?
            - run with nwalkers=40 but nstep  doubled to 400
            - run with nwalkers=40 but burnin doubled to 800

    - run-dg03r03 nwalkers=40, nstep=400 (doubled nstep, burnin same)
        - looks similar to r02
    - run-dg03r04 nwalkers=40, burnin=800 (doubled burnin, nstep same)
        - looks in between

    - so these are probably burnt in, but more sampling was needed. Real data
      may not burn in as quickly. options:
        - nwalkers=40 burnin=400 nstep=400 known bias
        - nwalkers=80 burnin=400 nstep=200 known bias
        - nwalkers=80 burnin=200 nstep=200 even less burnin
        - nwalkers=40 burnin=400 nstep=800 see if more nstep helps.
            - I think it is more important at this stage to look for better
              accuracy.  Let's do it

    - run-dg03r05 nwalkers=40, burnin=40, nstep=800
        - looks about the same


    - + run-dg03r06
        using prior after
        RUNNING

- sim-eg01 shear 0.08
    -  run-eg01r01 nwalkers=80, burnin=400,nstep=200
    -  run-eg01r02
        draw from priors, pretty high error, just a test
    -  run-eg01r03
        draw from priors, and more nstep

    - + run-eg01r04
        using prior after
        RUNNING

- sim-eg02 shear 0.04
    - trying to see if improves over 0.08, might be a clue
    - + run-eg02r01 nwalkers=80, burnin=400,nstep=200
        - looks quite similar to run-dg01r01

- sim-eg03
    - shear 0.08, size ratio sqrt(2)
    - run-eg03r01, nwalkers=80, rest standard
        - TERRIBLE.  Always low.  Maybe this is from a bad guess?  But that is
          hard to believe.
    - run-eg03r02 Draw guess from priors.
        - Better, but the exp are looking like shit.  What is going on?

- sim-dg04
    - shear 0.08, size ratio sqrt(2)
    - run-dg04r01, nwalkers=80, rest standard
        Looks OK but should try again with draw guess from priors

- going to run with shear=0.08 for now, nwalkers=80
    - we already have these runs
        - run-eg01r01
        - run-dg03r02
    - add these
        - run-dg03r06 - hold off on this
        - run-eg01r02 - hold off on this
        - run-eg03r01
        - run-dg03r04

- sim-gg01 shear=0.08
    - run-gg01r01 all standard 80,400,200
    - run-gg01r02 importance sampling
        - THIS ONE WAS BUGGY. DO NOT USE
    - run-gg01r03 
        - conjecture is that we need more samples not to get the tails but to
          see the spikiness of the prior.  And the mcmc sampler may not be able
          to sample it?  Who knows.
        - so try applying the prior after to see if we can get some of the
          spikiness
        - not remarkably better, but noiser and might need to sample more to
          see spikiness of prior.

    - run-gg01r04 increase nstep to 800 to see if I can get more spikiness and
      improve shear recovery.  Expect of order 4 hours to finish, 10pm.


- sim-gg04 shear=0.08, sigma ratio 1.0
    - run-gg04rtest3 not doing ring, just to see
        - not enough precision.  maybe come back to it
    - run-gg04r07 try expanding area over which objects are rendered and like
      calculated, etc.
         - rcas2259 slow again
    - also try larger render region.  Maybe "5-sigma" isn't general when the
      object is comparable to psf size?
          - run-gg07r01 not enough stats, but doesn't look good
          - run-gg07r02 more stats.  Same.

- sim-gg05 sigma ratio 1.0 but larger psf
    - same
- sim-gg06 sigma ratio 1.0 but even larger psf
    - same
- sim-gg07 larger render region (8 instead of 5 "sigma")
    - same
- sim-gg08 don't expand true shear; to test if this is because of numerical
  deriv.  This means low shear, 0.01

- near final run types
    - sigma ratio 2
        - gg01r01, but note *did* use priors during. also had error estimate
          wrong so it is actually much lower error.

        - gg01r04 more nstep and no g prior during
        - gg01r08
            - gg01rcomb

        - eg01r04
        - eg01r06
            - eg01rcomb

        - dg03r06
        - dg03r08*

    - sigma ratio 1.4
        - gg03r01
        - eg03r03
            - looks biased low
        - dg04r02
    - sigma ratio 1.0
        - sim-gg04 
            - gg04r01,gg04r02
                - was lower s/n than expected so ran two
                - looks biased.  Maybe error in psf fit affects these smaller
                  objects more?
                    - tried 2 gauss with gg04r03 and looks awful!
                - maybe psf not fully converged?  Use tighter tolerance and
                  higher maxiter?  Did a quick try and psf looks quite close.
                - maybe not burning in?
            - gg04r03
                - same parameters as gg04r01 but trying psf_ngauss=2
                - wow, that looks terrible.
            - gg04r04
                - doubling walkers to 160. Not better. deleted
                - draw truth.  Looks a bit better *maybe*
            - gg04r05
                - mca_a 2.0 didn't help
                - try mca_a 4.0, gives arate ~0.27.  Similar looking.
            - gg04r06
                - try same walkers but double burnin and nstep; doubtful
            - gg04rtest
                - nwalkers 400, no change!
                - try g prior during.  Same
        - sim-eg04
            - eg04r01,eg04r02.  Wierd oscillation in there.  I see it in the
              gg04 stuff as well.
            - eg04r03
        - sim-dg05
            - dg05r01
            - dg05r02
            - dg05r03 - running
            - dg05r04 - TODO

- pop idea
    - measurements
        - N galaxies
        - we have some set of measured parameters; this includes image noise.
    - assumptions
        - assume we know the PSF
        - assume we know the unlensed, intrinsic distribution of all the
          parameters for the underlying population.
        - assume that if a set of galaxies has the same distribution of
          parameters then it has the same shear.
    - procedure 

        - simulate a population of Ns galaxies drawn from the priors, and sheared
          by some specified amount.

        - Test this population of parameters against the one measured from
          the distribution of measured parameters (properly normalized).

        - Repeat and ask which shear value matches better.  Could use a
          non-linear fitter for this.  Take a guess from some simple shear
          estimate.

    - In principle, could use simple parameters like weighted moments.  Could
      even possibly use the *observed* moments rather than some fit?

    - for now probably test simple gaussian simulations, and do a standard fit
      and compare those.

    - the number of model evaluations would be Ns times the number of
      evaluations to get the estimator times the number steps needed to find
      the maximum likelihood.  For Ns=N, 100 steps in a non-linear for the
      estimator, and 100 steps to find the best shear, this would be say
      100*100*N, which is about the same as for the MCMC chain.  Adaptive
      moments would be fewer by probably a factor of ten.  Straight circular
      moments by a factor of 100.


    - practical difficulty is how we do the comparison?  Using a histogram in
      n-dimensions may not be practical.   Can we just look at the mean of some
      ellipticity parameters?

- look at universality of s/n bias vs T_s2n in a maximum likelihood fit or
  expectation value.
      - checked expectation value with priors and dev is quite far off the
        pattern.


- LM
    - run-gg01r05 
    - run-eg01r05
    - run-dg03r07
        - file too big

    - with randomization
    - gg01r06
        - only 1 random realization
    - gg01r07
        - 10 random realizations
        - errors at high s/n. Maybe lm pars?

    - need to get noise right; when calculating s2n_w I'm using all duplicate
      images.  And what about error for the combo over 10 images?  Currenly
      just multiplying by 0.5 inside code; better expand errors before sending
      I think, and let internal code be agnostic.

- metafit idea
    - Do a max like fit, or some other measure, to get a set of "observables".
    - Run a second fit (with guess from max like) to find the intrinsic object
      that would give those observables.  Each step here involves simulating
      the object convolved with the psf, rendering the image, and finding
      the observables.
    - problem  at low s/n the best fit will move all over the place.  This is
      probably prohibitive.
    - nice properties
        - should remove noise bias
        - fully incorporates all the observational details of the data such as
          noise, masked pixels, etc. in the simulation.
    - performance
        - typical number of steps for max like fit is about 50
        - thus a typical number of steps will be 50*50=2500
    - issues
        - failure of initial fit
        - poor errors on initial fit.  At very least we need to soften the
          covariance matrix.

    - gg01rtest
        - tried noise type "diff", difference of image and model.
            - total crap
        - noise type "random"
    - config/run-eg04rtest3.yaml
    - config/run-eg14rtest1.yaml
        - cosmos shapes
    - config/run-eg16rtest1.yaml
        -  nsub 1, why?

    - config/run-eg04mcal01
        - using true pars and nsub=16
        - s/n=1000
        - 0.66% error
        - trying again with different lm pars
            ftol: 1.0e-6
            xtol: 1.0e-3
            maxfev: 4000
        - Same
        - try no prior during
            - same
        - try guessing truth, random truth
            - same
        - try using fit instead of truth for metacal.  Just curious
            - blew up!  why?
        - run something stupid like em with no deconvolution?

    - config/run-eg04mcal02
        - s/n=25
        - using true pars and nsub=16
            - 1% error
    - issue: on real data we only have the pixelized psf so how
        does eric deal with that for convolution?

    - Anze points out this is estimator may be breaking down at high shear
    - sim-eg17
        - shear 0.01
        - big run: run-eg17mcal04.  frac error 0.002-0.004





# old shapesim stuff
- nsim-eg01
    - exp, sigma ratio 2 (T ratio 4)
    - ngmix-eg01r{01-12} nwalkers 40,burnin 400, nstep 200, 20 s2n bins
    - ngmix-eg01r{17-32} nwalkers 40,burnin 400, nstep 200, 11 s2n bins
    - ngmix-eg01r33 nwalkers 80,burnin 400, nstep 200, 11 s2n bins
        - with good s/n overall, but split. (no longer need to make
          multiple runs).  to test if more walkers helps....
          It looks pretty much the same as above.
- nsim-dg01
    - dev, sigma ratio 2 (T ratio 4)
    - ngmix-dg01r{01-16} nwalkers 40,burnin 400, nstep 200, 11 s2n bins
      fewer splits, so s/n in each is only ~1.0e-4

        - saw large outliers at first, but made guess tight ball around
          true and this went away. Not seen for exp....

    - ngmix-dg01r{17-32} more of the same

    - run-dg01r33
        - nwalkers=20 but new guess, to see if we get the same improvement.

- Things look much better since going to nwalkers=40 and taking a good
  guess.   exp that was already quite good before at this size ratio,
  so it is not much help in the interpretation.  But for dev I saw big
  outliers at intermediate s/n at first even with nwalkers=40, went
  away with better *guess*.
    - Interpretation
        - better guess helped with burnin, got rid of huge outliers for
          dev.  Did not change exp
        - better sampling after burnin generally improved result at low
          s/n
        - this means we might simulate the "better guess" by using
          more burnin.  To get better burnin *and* better sampling,
          would just increase the number of walkers.

    - tests of interpretation
        - use new style guess on dev but with 20 walkers.  If looks
          good then it might not be the sampling of the tails after
          burnin that was the fixer, just the burnin (or good guess).

