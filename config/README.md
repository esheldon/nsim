# vim: set filetype=markdown :

- lm
    - in tester in fitting.py I see no bias, what is different?
        - I wasn't simulating with sub-pixel integration
        - still doesn't explain difference between computers.
        - ran with nsub=1 in nsim and still seeing same thing
        - maybe fixed size box?

    - run-eg14r02 regular lm cosmos g dist T=2*Tpsf
        - bias not so bad, looks like im3shape

    - run-eg04rtest3 metacal ba dist T=Tpsf s/n=23 (27)
        - lots of negative bias, -0.027
        - either becuase smaller, becuase g dist is ba, or a bug with metacal version
        - the g_sens goes the wrong way! <g_sens> = 1.06!

        - re-running with draw-priors guess
            - same

    - run-eg04rtest4 regular lm ba g dist T=Tpsf
        - see if bug
            - it is very similar

    - run-eg15rtest1
        - same as run-eg14r02 but T=Tpsf
        - similar bias. So it is the size.

    - run-eg16rtest1
        - try nsub=1 : same

        - trying to match exactly: constant box size 25, T=4, nsub=1,
          shear=0.04, cen_sigma=0.5
          - now -2% bias!

- mh
    - eg12r02,eg12r03,eg12rtest6
        - working well except at s/n=10 for guess from maxlike (priors on all)
          and burnin=5000,nstep=5000
        - much faster obviously than the emcee
        - if I increase steps after burnin an improvement is ambiguous: it
          could effectively be acting like more burnin.
        - if improvement with more burning can safely leave nstep alone
        - if none it may be we can just increase the after burning steps
        - first try increaseing burnin to 10000.
    - eg12rtest7 10000,5000
        - no improvement
    - eg12rtest8 5000,10000
        - some improvement

    - might have had bad step sizes at low s/n due to bug in how the prio
        was used to set fdiff in lm
        - eg12rtest9 is same pars as eg12rtest6, see if it improves.
            - yes!!

        - still very small bias at s/n=10, ~8e-4 level
            - try temp
            - first s/n=100 with eg12rtest10

            - eg12r05 looks pretty good, although really just same bias on the
              positive side
    - eg04r04 - temp=1.5, Tobj=4 (fwhm ratio 1.2)
        - looks biased at s/n=10
    - eg04r05 - Tobj=4 temp=2.0

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

