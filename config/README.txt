REAL METACAL
--------------

important sims for paper
--------------------------

bd29mcal01
    - r50 match cosmos, psf fwhm match DES, ellip dist matching DES
bd29mcal02
    - same but with noise rotation; can combine with above

bd29mcal-t02, bd29mcal02

bdmask01mcal01
    - same as bd29 but with masking
bdmask01mcal02
    - same as bd29 but with masking

bdmask03mcal01
    - just bad pixels from DES
          m1: -4.615e-04 +/- 6.542e-04 c1: -1.720e-05 +/- 2.273e-05  r1: -0.00836
          m2: -7.305e-04 +/- 6.516e-04 c2: 7.453e-05 +/- 2.345e-05  r2: -0.033
          m:  -6.006e-04 +/- 4.620e-04
          c1: -1.716e-05 +/- 2.310e-05
          c2: 7.438e-05 +/- 2.311e-05
    - max like
          m1: -4.955e-01 +/- 2.887e-04 c1: 8.638e-06 +/- 1.003e-05  r1: -0.00836
          m2: -4.959e-01 +/- 2.814e-04 c2: 1.508e-03 +/- 1.013e-05  r2: -0.033
          m:  -4.957e-01 +/- 2.017e-04
          c1: 8.706e-06 +/- 1.009e-05
          c2: 1.507e-03 +/- 1.009e-05


bdmask04mcal01
    - just bad pixels from DES, with rotated added (not run)
    - not run

bd29stars-mcal01
    - 20,000,000 stars
bd29stars-mcal02
    - 20,000,000 stars, but using new code might have been broken

older
bd16: slightly undersampled (fwhm=0.8?), round psf. small galaxies
bd22: same as bd16 but non-round PSF
bd15: fwhm=0.9, larger galaxies, non-round psf
bd21: same as bd16 but non-round PSF and with masking

bd23: psf sizes drawn from examples from DES, gal sizes similar to COSMOS
bd24: same as 23 but gauss psf.
    - also added round size



new style simgs
=================

- specify distribution of flux not s2n
- noise specified in the run not the sim

- sim-ggnr03
- sim-ggnr03z
    - shear [0.08,0]
    - psf shape [0,0.05]
    - adjusted so that, for noise=1.0, mean s2n is about 12.5 with
        scatter 5

    - run-ggnr03zmcal-degrade01
        - noise 0.001
        - nrand 360

    - run-ggnr03mcal-t01
        - noise 1.0
        - 2,000,000
            meas: 0.0798909 +/- 0.000171054, -8.26837e-05 +/- 0.000168065
            fracdiff: -1.36e-03 +/- 2.14e-03

    - run-ggnr03mcal-01
        - noise 1.0
        - 20,000,000
            meas: 0.0794356 +/- 5.24952e-05, -0.000139687 +/- 5.5112e-05
            fracdiff: -7.05e-03 +/- 6.56e-04
        - this is 2.7 sigma from the previous one

    - run-ggnr03mcal-t01 revisit
        - run at slac and see if we get a bias
        - whoops, ran t01 instead
            - yep, much bigger bias. What happened?
    - run-ggnr03mcal-t02
        - at bnl

        meas: 0.0793634 +/- 0.00017042, -0.000232593 +/- 0.00017038
        fracdiff: -7.96e-03 +/- 2.13e-03

    - low s/n runs
        - noise 0.1 
        - run-ggnr03zmcal-degrade02
        - run-ggnr03mcal-t03
            meas: 0.0799988 +/- 0.000129255, -0.000120184 +/- 0.00011995
            fracdiff: -1.55e-05 +/- 1.62e-03

        - I was surprised that the mean came out fracdiff -0.0258 before
            correction, at <s/n>=125.  Maybe related to the prior bug?
            - yes I did a deep run and it was OK (sens=1)

    - had a bug in the prior, assumed linear transformation between
        r50 and T but there is a square
        - still, metacal should have fixed it.  Maybe the priors
            need to be broader
        - or maybe the "same noise" feature of the deep data only
            works if the priors are accurate
        - maybe will need to do not-same-noise deep data, with even
            more randoms. can try 160 for now but note it will be much slower,
            since we are limited by the metacal step ( for same noise we add noise
            to the metacal'd images).

    - run-ggnr03zmcal-degrade03
        - still noise=0.1
        - slac
    - run-ggnr03mcal-t04
        - still noise=0.1
        - bnl
        - the T distribution is wrong, pushed to higher values

    - run-ggnr03mcal-t05
        - noise=0.001
        - broad T prior
        - bnl

    - forgot that galsim takes the size of the *round* object!  This is
        why the measured values at high s/n are different than I would
        predict

    - run-ggnr03zmcal-degrade04
        - bnl (git not working at slac)
        - noise 1.0
        - broader T prior

    - run-ggnr03mcal-t06
        - using deep run-ggnr03zmcal-degrade04
        - bnl
        - noise 1.0
        - broader T prior
            meas: 0.0800078 +/- 0.00016898, -8.85498e-05 +/- 0.00017174
            fracdiff: 9.75e-05 +/- 2.11e-03
    - run-ggnr03mcal-t07
        - using deep run-ggnr03zmcal-degrade04
        - same as t06
        meas: 0.0798671 +/- 0.000159583, 0.000103818 +/- 0.000167474
        fracdiff: -1.66e-03 +/- 1.99e-03

    - run-ggnr03mcal-02
        - bnl
        - using deep run-ggnr03zmcal-degrade04
        - same as run-ggnr03mcal-t06 but longer run
            meas: 0.0797925 +/- 5.2574e-05, -2.88243e-05 +/- 5.28331e-05
            fracdiff: -2.59e-03 +/- 6.57e-04
          weighting according to mcal_s2n_r
            meas: 0.0797679 +/- 5.25903e-05, 2.58094e-05 +/- 5.28812e-05
            fracdiff: -2.90e-03 +/- 6.57e-04
        - averaging g and mcal_g
            fracdiff 6.09e-04 +/- 6.57e-04

    - run ggnr03zmcal-degrade05
        - bnl
        - prior ba width 0.3
    - run-ggnr03mcal-03
        - bnl
        - prior ba width 0.3

        - using deep ggnr03zmcal-degrade05
            meas: 0.0798205 +/- 5.29889e-05, 3.27795e-05 +/- 5.25502e-05
            fracdiff: -2.24e-03 +/- 6.62e-04
          weighting according to mcal_s2n_r

            meas: 0.0798439 +/- 5.29084e-05, 2.167e-05 +/- 5.24354e-05
            fracdiff: -1.95e-03 +/- 6.61e-04
        - averaging g and mcal_g
            fracdiff 4.46e-05 +/- 6.61e-04

        - using deep ggnr03zmcal-degrade06 with 0.1 noise

            meas: 0.0799327 +/- 5.30626e-05, -0.000135928 +/- 5.26465e-05
            fracdiff: -8.41e-04 +/- 6.63e-04

          interesting.  And I used a factor of ten like this in my
          other sims too, so maybe that helps.  But note the additive
          blew up!  Maybe additive is better measured starting with less
          noise?
        - with mcal_s2n_r weights

            meas: 0.0799583 +/- 5.29835e-05, -0.000145665 +/- 5.25334e-05
            fracdiff: -5.21e-04 +/- 6.62e-04

    - run-ggnr03mcal-04
        - same as run-ggnr03mcal-03 to get more stats

        - using deep ggnr03zmcal-degrade05
            - no weights
                meas: 0.0798354 +/- 5.36289e-05, 1.65783e-05 +/- 5.37265e-05
                fracdiff: -2.06e-03 +/- 6.70e-04

        - using deep ggnr03zmcal-degrade06 with 0.1 noise
            - no weights

            meas: 0.0799476 +/- 5.37044e-05, -0.000152168 +/- 5.38255e-05
            fracdiff: -6.55e-04 +/- 6.71e-04

            - mcal_s2n_r based weights
                meas: 0.0799773 +/- 5.3624e-05, -0.000162499 +/- 5.37394e-05
                fracdiff: -2.84e-04 +/- 6.70e-04

    - combine 03 and 04
        - weights
        fracdiff: -4.03e-04 +/- 4.71e-04

Now for a tougher test
But since gary used e not g for psf, I have set
the psf shape to [0,0.025]
- sim-bd06
- sim-bd06z
    - elliptical gauss psf
    - no shifts of center or bulge

    - run-bd06max-lownoise01    
        - exp
        - for getting T prior
        - used bad prior

    - high s/n runs
        x run-bd06zmcal-degrade01
            - don't use
        x run-bd06mcal-t01
            - don't use

        - run-bd06zmcal-degrade02
            - intentionally broad prior
        - run-bd06mcal-t02
            - intentionally broad prior
            - using run-bd06zmcal-degrade02
                meas: 0.0798944 +/- 0.000130532, 0.000160121 +/- 0.000123722
                fracdiff: -1.32e-03 +/- 1.63e-03
                weighting (but note noise not tuned for this sim/gal model etc)
                meas: 0.080089 +/- 0.000193841, 0.000230507 +/- 0.000192006
                fracdiff: 1.11e-03 +/- 2.42e-03

            - using run-bd06mcal-degrade01
                meas: 0.0803432 +/- 0.000131261, 0.000111218 +/- 0.000122873
                fracdiff: 4.29e-03 +/- 1.64e-03
            - using run-bd06mcal-degrade02
                meas: 0.0803638 +/- 0.000131301, 0.000180534 +/- 0.000122859
                fracdiff: 4.55e-03 +/- 1.64e-03

    - the response will be slightly different in the zero shear sim because it has
      different ellipticities (not sheared!) and different sizes (not sheared!)
        - we expect the sensitivity to be overestimated, since less ellip.
            in the deep sims.
        (- or I could have it completely backward: the reason we can't use the
            sensitivities from the regular data is because it is sheared)
        - might be somewhat overcome but using broad priors?  Not sure.
        - so the right thing to do in this unrealistic, constant shear situation
            is to shear the sim used to calculate the response
        - in a realistic sim with variable shear, we would also  have a deep field
            with random variable shear and all would work out
        - the whole noise adjustment thing was probably really making up
            for this fact
        - maybe shears of 0.08 but different random directions would be better?

    - run-bd06mcal-degrade01
        - this one using sheared sim, noise 0.001 target 0.01
            - doesn't work
    - run-bd06mcal-degrade02
        - this one using sheared sim, noise 0.0001 target 0.01
            - doesn't work

    - low s/n runs
        - run-bd06zmcal-degrade03
            - 0.1 target 1.0
        - run-bd06mcal-01
            - noise 1.0
            - using deep run-bd06zmcal-degrade03
                no weights
                    meas: 0.0799629 +/- 6.27301e-05, 0.000105573 +/- 6.25775e-05
                    fracdiff: -4.64e-04 +/- 7.84e-04
                old messed up weights
                    These weights don't even make sense
                    meas: 0.0799822 +/- 6.24348e-05, 9.36009e-05 +/- 6.22417e-05
                    fracdiff: -2.22e-04 +/- 7.80e-04
                new mcal_s2n_r weights
                    SN=0.18
                    meas: 0.0800931 +/- 6.25695e-05, 5.39724e-05 +/- 6.21891e-05
                    fracdiff: 1.16e-03 +/- 7.82e-04
                    SN=0.24 (a bit high)
                    meas: 0.0800773 +/- 6.23675e-05, 6.94577e-05 +/- 6.20595e-05
                    fracdiff: 9.66e-04 +/- 7.80e-04
        - run-bd06mcal-02
            - same as run-bd06mcal-01 for more stats
            - SLAC
                meas: 0.0798913 +/- 6.16718e-05, 8.2801e-05 +/- 6.24919e-05
                fracdiff: -1.36e-03 +/- 7.71e-04
        - combined run-bd06mcal-01, run-bd06mcal02
            fracdiff: -9.11e-04 +/- 5.50e-04
            with reasonable weighting
            fracdiff: 6.98e-04 +/- 5.48e-04
            with the crazy weighting
            fracdiff: -6.83e-04 +/- 5.47e-04


- sim-bd07
- sim-bd07z
    - elliptical moffat psf
    - no shifts of center or bulge
    - 3 gauss psf is considerably slower

    - high s/n
        - run-bd07zmcal-degrade01
            - noise 0.001, target 0.01
            * running at SLAC
        - run-bd07mcal-t01
            - relatively small sim
            - noise 0.01
            - ran at bnl, copied to slac 
            - using deep from run-bd07zmcal-degrade01
                meas: 0.0800917 +/- 0.000126313, 0.000522912 +/- 0.000124464
                fracdiff: 1.15e-03 +/- 1.58e-03
    - low s/n
        - run-bd07zmcal-degrade02
            - noise 0.1 target 1.0
        - run-bd07mcal-01
            -noise 1.0
            - using run-bd07zmcal-degrade02
                meas: 0.0798003 +/- 6.31235e-05, 0.000269427 +/- 6.43947e-05
                fracdiff: -2.50e-03 +/- 0.79e-03
            - noise bias undercorrected?  Maybe because the deep run has
                more uniform noise?

        - run-bd07max-lownoise01
            - to get T and flux priors

        - run-bd07zmcal-degrade03
            - wide priors from run-bd07max-lownoise01
            - noise 0.1 target 1.0
        - run-bd07mcal-t02
            - relatively short run for quick test of broad priors
                meas: 0.0796824 +/- 0.000207188, 0.000277947 +/- 0.000194573
                fracdiff: -3.97e-03 +/- 2.59e-03

        - run-bd07zmcal-degrade04
            - whiten same seed
            - slac
        - run-bd07mcal-t03
            - whiten same seed
            - bnl
                terrible

- sim-bd08
- sim-bd08z
    - reasonable flux/s2n distribution
    - elliptical moffat psf
    - no shifts of center or bulge
    - 3 gauss psf is considerably slower
    
    - run-bd08zmax-lownoise01
        - for priors
        - logT looks good with 10 gauss
        - logF that spike is not going to play well with the max like fitter.
            - maybe try the flat prior for flux?
    - run-bd08zmcal-degrade01
        - running slac
        - used two-sided-erf prior on counts
        - used slightly wider than lownoise prior on T
    - run-bd08mcal-t01
        - running bnl
        - used two-sided-erf prior on counts
        - used slightly wider than lownoise prior on T
    - run-bd08mcal-01
        - larger run but same settings
            meas: 0.0799852 +/- 5.15758e-05, 0.000318454 +/- 5.20944e-05
            fracdiff: -1.85e-04 +/- 6.45e-04
              ran different deep with starting noise 0.001
                  (run-bd08zmcal-degrade02)
                  meas: 0.0799569 +/- 5.15575e-05, 0.000382788 +/- 5.21092e-05
                  fracdiff: -5.39e-04 +/- 6.44e-04
              consistent, so that doesn't help leakage

    - run-bd08mcal-02
        - same as -01, looking for repeatability
            meas: 0.0800246 +/- 4.97376e-05, 0.000327837 +/- 5.18478e-05
            fracdiff: 3.07e-04 +/- 6.22e-04

    - 01 and02 combined
        meas: 0.0800049 +/- 3.58018e-05, 0.000323148 +/- 3.67486e-05
        fracdiff: 6.11e-05 +/- 4.48e-04
    - any calibration issues we see are apparently not from psf modeling
    - but the additive persists

    - run-bd08zmcal-degrade02
        - doing a run with less starting noise to see if leakage improves
        - leakage blew up when I went to moffat...
        - will include simple s/n

- sim-bd09
- sim-bd09z
    - same as sim-bd08* but added dev offset
    - run-bd09zmax-lownoise01
        - for priors
    - run-bd09zmcal-degrade01
    - run-bd09mcal-t01
        - bnl
            meas: 0.0798043 +/- 0.00016482, 0.000741121 +/- 0.000162534
            fracdiff: -2.45e-03 +/- 2.06e-03

            additive even higher
    - run-bd09mcal-01
        - slac 
            meas: 0.0800187 +/- 4.8071e-05, 0.000606897 +/- 4.93815e-05
            fracdiff: 2.34e-04 +/- 6.01e-04

        - try cuts or weights to help with leakage. Don't expect it
            to help because bright galaxies are no bigger!

            - cutting on noisy s2n_simple
                meas: 0.07887 +/- 5.09448e-05, 0.000597759 +/- 5.21425e-05
                fracdiff: -1.41e-02 +/- 6.37e-04
              note the deep did not get the same fraction removed, so
              something is fishy
            - weighting using s2n_simple
                meas: 0.0795862 +/- 5.17067e-05, 0.000672205 +/- 5.23488e-05
                fracdiff: -5.17e-03 +/- 6.46e-04

    - run-bd09zmcal-degrade02
        - using coellip4
    - run-bd09mcal-t02
        - using coellip4


- sim-bd10
- sim-bd10z
    - random orientation for psf
- run-bd10mcal-01
    - slac
        meas: 0.0800402 +/- 4.95431e-05, 0.000352045 +/- 4.76393e-05
        fracdiff: 5.03e-04 +/- 6.19e-04

    - bootstrapping responses
        meas: 0.0800415 +/- 5.21545e-05, 0.000362052 +/- 5.23894e-05
        fracdiff: 5.19e-04 +/- 6.52e-04

- sim-bd11
    - draw psf ellip from des-like distribution.  for des we see
        roughly a double gaussian
        psf g1 ~ -0.005 +/- 0.018
        psf g2 ~  0.007 +/- 0.018
      for this case I set the g1 offset to zero so we study primarily
      calibrations in g1 (shear is [0.08, 0.00]) and leakage in
      g2
    - run-bd11zmcal-degrade01
    - run-bd11mcal-01
        - straight averaging of signal, and using only the mean of
            sensitivities

            meas: 0.080009 +/- 4.86838e-05, 6.03224e-06 +/- 4.9374e-05
            fracdiff: 1.12e-04 +/- 6.09e-04

        - full bootstrap of both shallow and deep data simultaneously
            meas: 0.080009 +/- 5.03164e-05, 6.03224e-06 +/- 4.85065e-05
            fracdiff: 1.12e-04 +/- 6.29e-04

        - bootstrapping sensitivities (will be different for different bootstraps)

            meas: 0.0799963 +/- 5.03505e-05, -1.21533e-05 +/- 4.96312e-05
            fracdiff: -4.66e-05 +/- 6.29e-04
            meas: 0.0800196 +/- 5.02606e-05, 1.11246e-05 +/- 4.9608e-05
            fracdiff: 2.46e-04 +/- 6.28e-04
            meas: 0.0800049 +/- 5.04376e-05, -5.59096e-06 +/- 4.98377e-05
            fracdiff: 6.14e-05 +/- 6.30e-04
    - run-bd11mcal-02
        - not running with small update in ngmix to allow metacal obs sent
        - check repeatability and get more stats
            meas: 0.0800639 +/- 4.93544e-05, 0.000108922 +/- 4.9364e-05
            fracdiff: 7.99e-04 +/- 6.17e-04
        - full bootstrap of both shallow and deep data simultaneously
            meas: 0.0800639 +/- 5.04824e-05, 0.000108922 +/- 4.89934e-05
            fracdiff: 7.99e-04 +/- 6.31e-04

    - combined 01 and 02 (bootstrapped)
        meas: 0.0800526 +/- 4.10215e-05, 4.56466e-05 +/- 3.92758e-05
        fracdiff: 6.57e-04 +/- 5.13e-04

    - run-bd11mcal-03
        - as big as 01 and 02 combined
        - normal jackknifing
            meas: 0.0799976 +/- 3.47912e-05, 9.19536e-05 +/- 3.49262e-05
            fracdiff: -2.95e-05 +/- 4.35e-04
            !!!! turns out the error is really 2.3e-3 !!!!
            the bug in the dilate must have resulted in huge error in
            the psf corr.  This propagates into the 01 and 02 as well,
            and explains what seemed to be larger than calculated
            variance

        - full bootstrap
            meas: 0.0799976 +/- 3.65297e-05, 9.19536e-05 +/- 3.5082e-05
            fracdiff: -2.95e-05 +/- 4.57e-04

     - combined 01,02,03 (bootstrapped)
         meas: 0.0800172 +/- 2.55075e-05, 7.45221e-05 +/- 2.4586e-05
         fracdiff: 2.16e-04 +/- 3.19e-04

- run-bd11mcal-t02
    - fix in metacal, wasn't taking abs value of shear

- run-bd11mcal-t03
    - fix in metacal, don't interpolate twice

- run-bd11mcal-t05
    - currently for random tests

- run-bd11mcal-t06
    - new metacal fixes
    - ref run is bd12zmcal-degrade01


- sim-bd12
    - same as bd11 but with multiple shears
- run-bd12zmcal-degrade01
- run-bd12mcal-t01
    - using degrade01
        s1 m: -7.456e-03 +/- 7.441e-03 c: -2.801e-04 +/- 1.616e-04
        s2 m: -1.312e-02 +/- 6.685e-03 c: 4.028e-04 +/- 1.580e-04

        s1 m: -3.942e-03 +/- 7.145e-03 c: 5.976e-05 +/- 1.552e-04
        s2 m: 2.080e-03 +/- 6.406e-03 c: 2.326e-04 +/- 1.513e-04


    - using run-bd12zmcal-degrade02
        s1 m: -7.537e-04 +/- 7.161e-03 c: 6.191e-05 +/- 1.555e-04
        s2 m: 7.234e-03 +/- 6.435e-03 c: 2.116e-04 +/- 1.520e-04
    - using degrade02
        s1 m: -4.312e-03 +/- 7.463e-03 c: -2.788e-04 +/- 1.621e-04
        s2 m: -8.049e-03 +/- 6.718e-03 c: 3.818e-04 +/- 1.588e-04
    - using degrade03
        s1 m: -3.159e-03 +/- 7.141e-03 c: 6.245e-05 +/- 1.551e-04
        s2 m: 4.901e-03 +/- 6.422e-03 c: 2.230e-04 +/- 1.517e-04


- run-bd12mcal-t02
    - using reverted ngmix
    - using run-bd12zmcal-degrade01
        s1 m: -3.942e-03 +/- 7.145e-03 c: 5.976e-05 +/- 1.552e-04
        s2 m: 2.080e-03 +/- 6.406e-03 c: 2.326e-04 +/- 1.513e-0
    - using run-bd12zmcal-degrade02
        s1 m: -7.537e-04 +/- 7.161e-03 c: 6.191e-05 +/- 1.555e-04
        s2 m: 7.234e-03 +/- 6.435e-03 c: 2.116e-04 +/- 1.520e-04
    - using run-bd12zmcal-degrade03
        s1 m: -3.159e-03 +/- 7.141e-03 c: 6.245e-05 +/- 1.551e-04
        s2 m: 4.901e-03 +/- 6.422e-03 c: 2.230e-04 +/- 1.517e-04
    - using run-bd12zmcal-degrade04
        s1 m: -2.126e-03 +/- 7.159e-03 c: 5.961e-05 +/- 1.555e-04
        s2 m: 6.146e-03 +/- 6.430e-03 c: 2.230e-04 +/- 1.519e-04


- run-bd12zmcal-degrade02
    - revert ngmix (5efd6b394eece17738fdb99d870e100bc133b17c, "remove print")
- run-bd12zmcal-degrade03
    - ngmix 0649ad4ee7e319c44c768a4eb4b0e43868656df7
       "Merge pull request #23 from beckermr/master"
       did not clean the target dir before installing....

- run-bd12zmcal-degrade04
    - ngmix 858916e16359e89291a87d4da70dbd15b97fe442
        "fix bugs in symmetrize"

- run-bd12mcal-01
    - using degrade01
        s1 m: -5.155e-03 +/- 1.529e-03 c: 2.592e-05 +/- 3.321e-05
        s2 m: -5.927e-03 +/- 1.371e-03 c: 1.112e-04 +/- 3.243e-05
    - using degrade02
        s1 m: -1.964e-03 +/- 1.530e-03 c: 2.776e-05 +/- 3.323e-05
        s2 m: -7.927e-04 +/- 1.379e-03 c: 8.947e-05 +/- 3.260e-05
    - using degrade03
        s1 m: -4.364e-03 +/- 1.526e-03 c: 2.830e-05 +/- 3.316e-05
        s2 m: -3.118e-03 +/- 1.373e-03 c: 1.013e-04 +/- 3.247e-05
    - using degrade04
        s1 m: -3.342e-03 +/- 1.533e-03 c: 2.574e-05 +/- 3.329e-05
        s2 m: -1.887e-03 +/- 1.375e-03 c: 1.011e-04 +/- 3.252e-05

- run-bd12mcal-02
    - reverted ngmix "remove print" 5efd6b394eece17738fdb99d870e100bc133b17c
    - should go with run-bd12zmcal-degrade02
        s1 m: -2.693e-03 +/- 1.705e-03 c: -4.615e-05 +/- 3.703e-05
        s2 m: -2.562e-03 +/- 1.442e-03 c: -1.704e-05 +/- 3.410e-05

- run-bd12mcal-t03
    - using new Rnoise corrections
    - 2 million
    - using Rnoise
      m1: 8.532e-03 +/- 7.185e-03 c1: 9.459e-05 +/- 1.561e-04  r1: 0.0404
      m2: 3.866e-03 +/- 6.752e-03 c2: -4.637e-04 +/- 1.596e-04  r2: -0.0438
      Q: 439

    - using Rpsf_noise
      m1: 9.164e-03 +/- 7.170e-03 c1: 9.248e-05 +/- 1.558e-04  r1: 0.0404
      m2: 4.604e-03 +/- 6.719e-03 c2: -4.719e-04 +/- 1.589e-04  r2: -0.0438
      Q: 410

    - using mean of two
      m1: 8.899e-03 +/- 7.173e-03 c1: 9.433e-05 +/- 1.558e-04  r1: 0.0404
      m2: 4.305e-03 +/- 6.724e-03 c2: -4.685e-04 +/- 1.590e-04  r2: -0.0438
      Q: 421.692



- run-bd12mcal-04
    - using new Rnoise corrections
    - ? 40 million ?

- Using new T prior and reverted ngmix 
    - ngmix "remove print" 5efd6b394eece17738fdb99d870e100bc133b17c
    - T prior from run-bd12zmax-lownoise01

    - run-bd12zmcal-degrade05

    - run-bd12mcal-03
        - T prior from run-bd12zmax-lownoise01
        - run-bd12zmcal-degrade05
            s1 m: -4.767e-03 +/- 1.656e-03 c: -2.602e-05 +/- 3.598e-05
            s2 m: -4.349e-03 +/- 1.495e-03 c: 1.440e-04 +/- 3.535e-05

    - run-bd11mcal-04
        - using run-bd12zmcal-degrade05


tracking bug
--------------
ngmix/nsim branches btest. Basically master but stripped
down bootstrapper

- run-bd12zmcal-degrade06
- run-bd11mcal-t07
    meas: 0.0797489 +/- 0.000177108, -3.47196e-05 +/- 0.000148962
    s1 fracdiff: -3.14e-03 +/- 2.21e-03
- run-bd11mcal-t08
    - same config
    s1 fracdiff: -1.08e-03 +/- 2.10e-03


- run-bd12zmax-lownoise02
- run-bd12zmcal-degrade07
- run-bd11mcal-t09
    - small 2 million run
    meas: 0.0798965 +/- 0.000168158, 0.000196409 +/- 0.000147272
    s1 fracdiff: -1.29e-03 +/- 2.10e-03
- run-bd11mcal-05
    - 40 million
        meas: 0.0798987 +/- 4.63937e-05, 0.000109258 +/- 3.36177e-05
        s1 fracdiff: -1.27e-03 +/- 5.80e-04


- run-ggnr03mn-t01
    - first metanoise run
    - noise 0.1
    - nrand 100 -> noise_target sqrt(100)=10 times original noise
    - expect this is ten times more noisy than the degrade version
    - frac error about 0.007 
- run-ggnr03mn-t04
    - 10 times more
    - intended noise was 1.0
        whoops, ran with 0.1. Biased anyway

moments
- run-ggnr03mn-t02
    - moments, with metanoise
- run-ggnr03mn-t03
    - moments, just very high s/n but no metanoise
    - ran this first
        - horrible

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
        - bnl
        meas: 0.0352182 +/- 0.000116463, 0.000519263 +/- 0.000115394
        fracdiff: 6.24e-03 +/- 3.33e-03

        some additive is back

        sensitivity is much broader than bd02-04.  Either due to
        psf or due to centroid shift

        making cut s2n_r > 10 gives same answer!  must be a bug

    - run-bd01mcal-01
        - planned full run


- sim-bd13
    - round gaussian psf
    - run-bd13mcal-dt-t04
        - 2 million
        - detrending
           m1: 2.690e-04 +/- 4.417e-03 c1: -1.162e-04 +/- 1.589e-04  r1: -0.13
           m2: -1.758e-03 +/- 4.437e-03 c2: -1.511e-04 +/- 1.539e-04  r2: -0.00774
        - not removing offset
           m1: -2.117e-03 +/- 4.407e-03 c1: -1.160e-04 +/- 1.585e-04  r1: -0.13
           m2: -4.308e-03 +/- 4.426e-03 c2: -1.509e-04 +/- 1.536e-04  r2: -0.00774
    - run-bd13mcal-dt01 8 million
        Not subtracting offset
          m1: -8.332e-04 +/- 2.118e-03 c1: -7.281e-05 +/- 7.619e-05  r1: -0.13
          m2: -1.188e-03 +/- 2.305e-03 c2: 4.958e-06 +/- 7.998e-05  r2: -0.00774
        Subtracting offset
           m1: 1.653e-03 +/- 2.123e-03 c1: -7.295e-05 +/- 7.638e-05  r1: -0.13
           m2: 1.401e-03 +/- 2.311e-03 c2: 4.897e-06 +/- 8.019e-05  r2: -0.00774

    - short runs using less added noise, to see if we go through
        zero in delta R
        - we do, but the answer is biased using the A
        - run-bd13mcal-dt-t05
            m1: -1.594e-02 +/- 8.500e-03 c1: -1.774e-04 +/- 3.058e-04  r1: -0.13
            m2: -3.084e-02 +/- 8.288e-03 c2: 3.901e-04 +/- 2.876e-04  r2: -0.00774
        - run-bd13mcal-dt-t06
            m1: -3.196e-02 +/- 8.459e-03 c1: 2.129e-04 +/- 3.043e-04  r1: -0.13
            m2: -2.152e-02 +/- 8.790e-03 c2: 3.652e-04 +/- 3.050e-04  r2: -0.00774

    - big run, 100 million
            m1: 9.878e-04 +/- 5.780e-04 c1: 1.953e-06 +/- 2.079e-05  r1: -0.13
            m2: 1.644e-04 +/- 6.520e-04 c2: 7.446e-06 +/- 2.262e-05  r2: -0.00774
        some cuts
            m1:  6.087e-04 +/- 5.745e-04 c1: 1.052e-06 +/- 2.067e-05  r1: -0.13
            m2: -1.038e-04 +/- 6.531e-04 c2: 6.450e-06 +/- 2.266e-05  r2: -0.00774
        weights 1/(1+(10/s2n)^2 )
            m1: -1.572e-04 +/- 5.826e-04 c1: -5.324e-06 +/- 2.096e-05  r1: -0.13
            m2: -8.463e-04 +/- 6.578e-04 c2: 5.631e-06 +/- 2.282e-05  r2: -0.00774
        weights 1/(1+(20/s2n)^2 )
            m1: -5.831e-04 +/- 6.301e-04 c1: -1.418e-05 +/- 2.267e-05  r1: -0.13
            m2: -1.207e-03 +/- 6.927e-04 c2: 7.740e-06 +/- 2.404e-05  r2: -0.00774


- sim-bd15
    - same as bd12 but with everything a bit better resolved
    * run-bd15mcal-dt01
      half way
      Not subtracting offset
        m1: -2.190e-03 +/- 2.387e-03 c1: -7.640e-05 +/- 8.586e-05  r1: -0.13
        m2: -2.544e-03 +/- 2.264e-03 c2: 2.559e-04 +/- 7.855e-05  r2: -0.00774
      Subtracting offset
        m1: 5.643e-04 +/- 2.393e-03 c1: -7.663e-05 +/- 8.610e-05  r1: -0.13
        m2: 1.016e-04 +/- 2.270e-03 c2: 2.386e-04 +/- 7.876e-05  r2: -0.00774

    - run-bd15mcal-dt-t01
        Purposefully using just a gaussian for psf
        Not subtracting offset
            m1: -5.206e-03 +/- 8.402e-03 c1: -4.138e-05 +/- 3.023e-04  r1: -0.13
            m2: -9.162e-03 +/- 9.199e-03 c2: 2.818e-04 +/- 3.192e-04  r2: -0.00774
        Subtracting offset
            m1: -3.017e-03 +/- 8.421e-03 c1: -4.143e-05 +/- 3.029e-04  r1: -0.13
            m2: -6.523e-03 +/- 9.224e-03 c2: 2.651e-04 +/- 3.200e-04  r2: -0.00774

    - adjusted noise to 0.76 to get mode of s/n pdf at 10
    - run-bd15mcal-dt02
        - gauss psf and gal
        m1: -4.603e-04 +/- 6.285e-04 c1: -2.080e-05 +/- 2.261e-05
        m2: -1.711e-03 +/- 6.455e-04 c2: 4.171e-05 +/- 2.240e-05
        m:  -1.069e-03 +/- 4.522e-04
        c1: -1.795e-05 +/- 2.260e-05
        c2: 4.154e-05 +/- 2.250e-05
    - run-bd15mcal-dt03
        - another new seed, otherwise the same as dt02

        m1: 4.106e-04 +/- 6.364e-04 c1: -1.432e-05 +/- 2.290e-05
        m2: 7.672e-04 +/- 6.640e-04 c2: -1.337e-05 +/- 2.304e-05
        m:  5.840e-04 +/- 4.608e-04
        c1: -1.513e-05 +/- 2.303e-05
        c2: -1.333e-05 +/- 2.293e-05

    - combined dt02,dt03
        m1: -5.198e-05 +/- 4.386e-04 c1: -1.595e-05 +/- 1.578e-05
        m2: -5.917e-04 +/- 4.684e-04 c2: 1.775e-05 +/- 1.625e-05
        m:  -3.144e-04 +/- 3.216e-04
        c1: -1.472e-05 +/- 1.607e-05
        c2: 1.768e-05 +/- 1.600e-05

        m1: -1.847e-05 +/- 4.412e-04 c1: -1.749e-05 +/- 1.587e-05
        m2: -4.566e-04 +/- 4.583e-04 c2: 1.382e-05 +/- 1.590e-05
        m:  -2.315e-04 +/- 3.189e-04
        c1: -1.649e-05 +/- 1.594e-05
        c2: 1.376e-05 +/- 1.587e-05


- sim-bd16
    - same as bd12 but psf is round on average
    - run-bd16mcal-dt-t01
        - 500,000
        - looks fine
    - run-bd16mcal-dt-t02
        - messed up the noise: this is high s/n
    - run-bd16mcal-dt01
        - 100 million
        - noise tuned for mode closer to 10
        - final
            m1: 6.089e-04 +/- 6.087e-04 c1: -3.087e-05 +/- 2.145e-05  r1: 0.00386
            m2: 3.114e-04 +/- 6.442e-04 c2: 2.152e-05 +/- 2.128e-05  r2: -0.0233

            m: 4.685e-04 +/- 4.434e-04
            c1: -3.093e-05 +/- 2.143e-05
            c2: 2.146e-05 +/- 2.141e-05

    - run-bd16mcal-dt01 + run-bd18mcal-dt01
        - adding 10% stars
            m1: 6.057e-04 +/- 6.512e-04 c1: -2.360e-05 +/- 2.295e-05  r1: 0.00386
            m2: 3.434e-04 +/- 6.681e-04 c2: 2.356e-05 +/- 2.207e-05  r2: -0.0233

            m: 4.813e-04 +/- 4.679e-04
            c1: -2.367e-05 +/- 2.258e-05
            c2: 2.361e-05 +/- 2.256e-05

    - selections
         note using different fake_shear values (0.01 vs 0.032 or whatever)
        did not make a difference.  Using the mean R etc. didn't work
        at all.

        running now without subtracting mean, but I think this didn't work
        better either.  also it is needed to remove net psf anisotropy
        anyway...
        - looks better for 10

        what else to try to get this under 1.0e-3?

        - -1.0 < R < 1.2 and no correction for selection
        removes 0.4%
        m1: -4.180e-03 +/- 6.022e-04 c1: -3.062e-05 +/- 2.122e-05
        m2: -4.323e-03 +/- 6.415e-04 c2: 2.184e-05 +/- 2.119e-05
        m:  -4.247e-03 +/- 4.401e-04
        c1: -3.063e-05 +/- 2.126e-05
        c2: 2.178e-05 +/- 2.126e-05


        - s/n > 20
            - no correction
                m: -2.662e-03 +/- 5.940e-04
                c1: -2.766e-05 +/- 2.869e-05
                c2: -2.223e-05 +/- 2.869e-05
            - with correction
                 m: -8.741e-04 +/- 5.947e-04
                 c1: -2.775e-05 +/- 2.873e-05
                 c2: -2.247e-05 +/- 2.873e-05
                 # not subtracting mean m = 1.142e-03 +/- 6.048e-04
        - s/n > 15 
            - no correction
                m: -2.315e-03 +/- 5.324e-04
                c1: -4.342e-06 +/- 2.571e-05
                c2: 2.091e-06 +/- 2.572e-05
            - with correction
                m: -9.060e-04 +/- 5.328e-04
                c1: -4.375e-06 +/- 2.574e-05
                c2: 1.928e-06 +/- 2.574e-05
        - s/n > 10
            - no correction
                m: -1.520e-03 +/- 4.616e-04
                c1: -2.713e-05 +/- 2.230e-05
                c2: 1.718e-05 +/- 2.230e-05
            - with correction
                m: -1.040e-03 +/- 4.618e-04
                c1: -2.716e-05 +/- 2.231e-05
                c2: 1.709e-05 +/- 2.231e-05
            - new correction

        - s/n > 7.5
            - with correction
                m: -5.480e-04 +/- 4.448e-04
                c1: -3.217e-05 +/- 2.149e-05
                c2: 2.183e-05 +/- 2.149e-05
        - s/n > 5
            - no correction
                m: 3.047e-04 +/- 4.436e-04
                c1: -3.081e-05 +/- 2.143e-05
                c2: 2.103e-05 +/- 2.143e-05
            - with correction
                m: 3.335e-04 +/- 4.436e-04
                c1: -3.082e-05 +/- 2.143e-05
                c2: 2.100e-05 +/- 2.143e-05

    - simple s2n weights
        - no corr
        s2n_soft 10.0
        m1: -1.009e-03 +/- 6.067e-04 c1: -2.622e-05 +/- 2.138e-05
        m2: -7.225e-04 +/- 6.289e-04 c2: 1.041e-05 +/- 2.077e-05
        m: -8.749e-04 +/- 4.375e-04
        c1: -2.620e-05 +/- 2.113e-05
        c2: 1.053e-05 +/- 2.114e-05

        - with correction
        m1: -2.367e-04 +/- 6.071e-04 c1: -2.624e-05 +/- 2.140e-05
        m2: -1.227e-03 +/- 6.286e-04 c2: 1.041e-05 +/- 2.076e-05
        m: -6.995e-04 +/- 4.380e-04
        c1: -2.630e-05 +/- 2.116e-05
        c2: 9.999e-06 +/- 2.116e-05

    - erf s2n weights 5

    - noise weights shapenoise=0.20
          m1: -2.757e-03 +/- 5.977e-04 c1: -2.587e-05 +/- 2.107e-05  r1: 0.00386
          m2: -2.901e-03 +/- 6.250e-04 c2: 1.595e-05 +/- 2.065e-05  r2: -0.0233
          m:  -2.824e-03 +/- 4.329e-04
          c1: -2.588e-05 +/- 2.091e-05
          c2: 1.589e-05 +/- 2.091e-05


- sim-bd17
    - same as bd12 but big, fixed psf anisitropy

- sim-bd18
    - stars
    - run-bd18mcal-dt-t01     500,000
    - run-bd18mcal-dt-t02     500,000 using flat priors and not log
    - run-bd18mcal-dt-t03   2,000,000
    - run-bd18mcal-dt01    10,000,000

- sim-bd19
    - same as bd17 but with a random bad pixel, and using noise
        to replace the value
    - run-bd19mcal-dt-t01 500,000
        m1: -3.231e-03 +/- 9.709e-03 c1: 2.089e-04 +/- 3.235e-04  r1: -0.00487
        m2: 2.193e-03 +/- 9.419e-03 c2: 3.631e-04 +/- 3.377e-04  r2: 0.103

sim-bd20
    - run-bd20mcal-dt-t01
        - added full masking, replaced by random noise
            - shitty
              m1: -4.522e-02 +/- 8.946e-03 c1: -1.930e-04 +/- 2.981e-04
              m2: -3.343e-02 +/- 8.518e-03 c2: 2.702e-04 +/- 3.055e-04
sim-bd21
    - significant psf anisotropy
    - run-bd21mcal-dt-t01
        - replaced by best fit, which is gaussian, so this might suck
              m1: 1.499e-02 +/- 9.326e-03 c1: -1.646e-04 +/- 3.107e-04
              m2: 4.376e-04 +/- 8.974e-03 c2: 4.824e-04 +/- 3.218e-04
    - run-bd21mcal-dt01
        - masked replaced by best fit model
        - 100 million run
             m1: -5.829e-04 +/- 7.125e-04 c1: -2.114e-06 +/- 2.374e-05
             m2: -9.088e-04 +/- 6.610e-04 c2: 1.152e-04 +/- 2.370e-05

             m: -7.574e-04 +/- 4.859e-04
             c1: -2.177e-06 +/- 2.372e-05
             c2: 1.156e-04 +/- 2.379e-05

        - no correction
             m1: 1.046e-01 +/- 7.871e-04 c1: -2.204e-06 +/- 2.623e-05
             m2: 1.044e-01 +/- 7.313e-04 c2: -5.596e-04 +/- 2.622e-05
             m: 1.045e-01 +/- 5.371e-04
             c1: -2.282e-06 +/- 2.622e-05
             c2: -5.593e-04 +/- 2.630e-05

        - select (-1.0 < R < 1.2) no correction
            m1: -5.217e-03 +/- 7.104e-04 c1: -2.803e-06 +/- 2.367e-05
            m2: -5.428e-03 +/- 6.594e-04 c2: 1.375e-04 +/- 2.365e-05
            m:  -5.330e-03 +/- 4.846e-04
            c1: -2.785e-06 +/- 2.365e-05
            c2: 1.379e-04 +/- 2.372e-05

        - select s2n_r > 10
              m1: -2.157e-03 +/- 7.307e-04 c1: 1.721e-05 +/- 2.434e-05  r1: -0.00487
              m2: -1.830e-03 +/- 7.148e-04 c2: 1.028e-04 +/- 2.563e-05  r2: 0.103
              m:  -1.982e-03 +/- 5.120e-04
              c1: 1.718e-05 +/- 2.499e-05
              c2: 1.022e-04 +/- 2.506e-05

        - s2n erf weighting with corr, rolloff at 10 width 5
            --reflect
              m1: -1.608e-03 +/- 7.024e-04 c1: 1.900e-05 +/- 2.340e-05
              m2: -1.605e-03 +/- 6.681e-04 c2: 1.139e-04 +/- 2.396e-05
              m:  -1.607e-03 +/- 4.851e-04
              c1: 1.900e-05 +/- 2.368e-05
              c2: 1.139e-04 +/- 2.375e-05

    - adding run-stars01mcal-dt01
        - corrected
            m1: -4.484e-04 +/- 7.235e-04 c1: 8.026e-06 +/- 2.411e-05
            m2: -9.286e-04 +/- 6.910e-04 c2: 1.302e-04 +/- 2.478e-05

            m: -7.048e-04 +/- 5.008e-04
            c1: 8.067e-06 +/- 2.444e-05
            c2: 1.310e-04 +/- 2.451e-05

            - no corr
                m1: 1.191e-01 +/- 8.099e-04 c1: 9.282e-06 +/- 2.699e-05
                m2: 1.187e-01 +/- 7.745e-04 c2: -6.364e-04 +/- 2.777e-05
                m: 1.189e-01 +/- 5.609e-04
                c1: 9.315e-06 +/- 2.738e-05
                c2: -6.358e-04 +/- 2.746e-05


    - run-bd21mcal-dt02
        - 2 gaussian psf model
            m1: -7.176e-04 +/- 7.025e-04 c1: -4.990e-05 +/- 2.341e-05
            m2: 7.741e-06 +/- 6.679e-04 c2: 9.373e-05 +/- 2.395e-05
            m: -3.302e-04 +/- 4.853e-04
            c1: -4.997e-05 +/- 2.368e-05
            c2: 9.248e-05 +/- 2.375e-05
        - no corr
            m1: 1.089e-01 +/- 7.799e-04 c1: -5.524e-05 +/- 2.599e-05
            m2: 1.095e-01 +/- 7.406e-04 c2: -6.111e-04 +/- 2.656e-05

            m: 1.092e-01 +/- 5.383e-04
            c1: -5.530e-05 +/- 2.627e-05
            c2: -6.121e-04 +/- 2.635e-05

        - s/n > 10
            m1: -1.432e-03 +/- 7.298e-04 c1: -6.217e-05 +/- 2.431e-05
            m2: -8.157e-04 +/- 6.777e-04 c2: 7.212e-05 +/- 2.430e-05
            m:  -1.103e-03 +/- 4.981e-04
            c1: -6.222e-05 +/- 2.431e-05
            c2: 7.105e-05 +/- 2.438e-05

        - s2n erf weighting with corr
            --reflect
            sel: [ 1.00051627  1.00050059]
            m1: -1.186e-03 +/- 7.076e-04 c1: -6.124e-05 +/- 2.358e-05
            m2: -5.064e-04 +/- 6.534e-04 c2: 7.602e-05 +/- 2.343e-05
            m: -8.229e-04 +/- 4.817e-04
            c1: -6.130e-05 +/- 2.351e-05
            c2: 7.485e-05 +/- 2.358e-05

    - run-bd21mcal-dt03
        - more stats to augment dt01
          m1: -5.298e-06 +/- 6.473e-04 c1: 8.556e-06 +/- 2.157e-05  r1: -0.00487
          m2: -4.141e-04 +/- 6.097e-04 c2: 7.335e-05 +/- 2.186e-05  r2: 0.103
          m:  -2.236e-04 +/- 4.449e-04
          c1: 8.591e-06 +/- 2.172e-05
          c2: 7.406e-05 +/- 2.178e-05

    - combined dt01 and dt03
          m1: -2.968e-04 +/- 4.979e-04 c1: 3.175e-06 +/- 1.659e-05  r1: -0.00487
          m2: -6.618e-04 +/- 4.625e-04 c2: 9.440e-05 +/- 1.658e-05  r2: 0.103
          m:  -4.917e-04 +/- 3.398e-04
          c1: 3.206e-06 +/- 1.659e-05
          c2: 9.503e-05 +/- 1.664e-05

sim-bd22
    - like bd21, significant psf anisotropy but no masking
    - run-bd22mcal-dt01
        - 100 million
            m1: -4.448e-05 +/- 6.640e-04 c1: 1.360e-05 +/- 2.340e-05
            m2: -6.225e-04 +/- 7.499e-04 c2: 1.103e-04 +/- 2.477e-05
            m:  -3.147e-04 +/- 5.002e-04
            c1: 1.356e-05 +/- 2.416e-05
            c2: 1.101e-04 +/- 2.416e-05
    - run-bd22mcal-dt02
        - 100 million
        - combined a bit early (1240 missing)
            (also some 430 didn't finish, running now)

    - combined run
        m1: 4.905e-04 +/- 4.531e-04 c1: 1.823e-05 +/- 1.597e-05
        m2: -9.254e-04 +/- 5.136e-04 c2: 1.135e-04 +/- 1.697e-05
        m:  -1.715e-04 +/- 3.432e-04
        c1: 1.814e-05 +/- 1.658e-05
        c2: 1.129e-04 +/- 1.658e-05

- sim-bd23
    - psf of psf fwhm drawn from DES
    - galaxy sizes similar to cosmos
    - run-bd23mcal-dt01
        m1: 4.288e-03 +/- 8.256e-04 c1: 2.904e-06 +/- 2.970e-05
        m2: 2.633e-03 +/- 7.849e-04 c2: 2.242e-05 +/- 2.723e-05
        m:  3.483e-03 +/- 5.724e-04
        c1: 6.670e-06 +/- 2.861e-05
        c2: 2.219e-05 +/- 2.848e-05
    - run-bd23mcal-dt02
        - using closer to "true" priors
        - 95/100 done
            m1: 2.217e-03 +/- 8.228e-04 c1: 4.833e-05 +/- 2.960e-05  r1: -0.13
            m2: 3.263e-03 +/- 7.853e-04 c2: 6.479e-05 +/- 2.725e-05  r2: -0.00774
            m:  2.725e-03 +/- 5.709e-04
            c1: 4.594e-05 +/- 2.853e-05
            c2: 6.494e-05 +/- 2.841e-05

        - 24/100 done
         m1: 2.161e-03 +/- 1.524e-03 c1: 2.247e-05 +/- 5.484e-05  r1: -0.13
         m2: 4.193e-03 +/- 1.539e-03 c2: 8.242e-05 +/- 5.340e-05  r2: -0.00774
         m:  3.149e-03 +/- 1.086e-03
         c1: 1.784e-05 +/- 5.430e-05
         c2: 8.270e-05 +/- 5.406e-05
        - 41/100 done
          m1: 2.818e-03 +/- 1.267e-03 c1: 5.731e-05 +/- 4.557e-05
          m2: 2.686e-03 +/- 1.170e-03 c2: 6.510e-05 +/- 4.061e-05
          m:  2.754e-03 +/- 8.654e-04
          c1: 5.761e-05 +/- 4.325e-05
          c2: 6.509e-05 +/- 4.306e-05

         after 41/100 -1.0 < R < 1.2 and no correction for selection
         removes 0.08%
          m1: 1.516e-03 +/- 1.265e-03 c1: 5.601e-05 +/- 4.549e-05  r1: -0.13
          m2: 1.126e-03 +/- 1.165e-03 c2: 6.833e-05 +/- 4.041e-05  r2: -0.00774
          m:  1.326e-03 +/- 8.627e-04
          c1: 5.690e-05 +/- 4.312e-05
          c2: 6.827e-05 +/- 4.293e-05

        - after 24/100

            s/n > 7.5
          m1: -4.657e-03 +/- 1.534e-03 c1: 2.454e-05 +/- 5.518e-05  r1: -0.13
          m2: -8.064e-04 +/- 1.524e-03 c2: 6.043e-05 +/- 5.286e-05  r2: -0.00774
          m:  -2.785e-03 +/- 1.087e-03
          c1: 1.578e-05 +/- 5.431e-05
          c2: 6.096e-05 +/- 5.407e-05

            s/n > 10
          m1: -5.250e-03 +/- 1.565e-03 c1: 4.900e-06 +/- 5.632e-05  r1: -0.13
          m2: -4.795e-05 +/- 1.607e-03 c2: 9.393e-05 +/- 5.576e-05  r2: -0.00774
          m:  -2.721e-03 +/- 1.129e-03
          c1: -6.941e-06 +/- 5.644e-05
          c2: 9.465e-05 +/- 5.619e-05

            s/n > 15 
          m1: -4.483e-03 +/- 1.577e-03 c1: -1.986e-05 +/- 5.672e-05
          m2: 1.164e-03 +/- 1.655e-03 c2: 8.841e-05 +/- 5.742e-05
          m:  -1.738e-03 +/- 1.151e-03
          c1: -3.272e-05 +/- 5.751e-05
          c2: 8.919e-05 +/- 5.726e-05
            s/n > 20
          m1: -8.778e-03 +/- 1.736e-03 c1: -2.244e-05 +/- 6.247e-05  r1: -0.13
          m2: -3.416e-03 +/- 1.789e-03 c2: 1.117e-04 +/- 6.208e-05  r2: -0.00774
          m:  -6.171e-03 +/- 1.254e-03
          c1: -3.465e-05 +/- 6.268e-05
          c2: 1.125e-04 +/- 6.240e-05

    - run-bd23mcal-dt-t05
        - closer to "true" priors
        - round T
    - run-bd23mcal-dt03
        - closer to "true" priors
        - round T

      m1: 4.449e-03 +/- 7.539e-04 c1: 2.415e-05 +/- 2.712e-05  r1: -0.13
      m2: 2.875e-03 +/- 8.049e-04 c2: 3.850e-05 +/- 2.793e-05  r2: -0.00774
      m:  3.684e-03 +/- 5.532e-04
      c1: 2.773e-05 +/- 2.765e-05
      c2: 3.828e-05 +/- 2.753e-05

     s/n > 10
      m1: -6.351e-03 +/- 7.639e-04 c1: 1.921e-05 +/- 2.748e-05  r1: -0.13
      m2: -7.442e-03 +/- 8.204e-04 c2: 1.071e-05 +/- 2.847e-05  r2: -0.00774
      m:  -6.882e-03 +/- 5.618e-04
      c1: 2.169e-05 +/- 2.808e-05
      c2: 1.056e-05 +/- 2.796e-05
    s/n > 20
      m1: -1.230e-02 +/- 8.517e-04 c1: 4.450e-06 +/- 3.064e-05  r1: -0.13
      m2: -1.250e-02 +/- 1.008e-03 c2: 1.694e-06 +/- 3.497e-05  r2: -0.00774
      m:  -1.240e-02 +/- 6.599e-04
      c1: 4.906e-06 +/- 3.298e-05
      c2: 1.666e-06 +/- 3.284e-05

- sim-bd24
    - gaussian psf, just to see if bias decreases relative to bd23
    - run-bd24mcal-dt01
        m1: 3.063e-03 +/- 7.314e-04 c1: 2.836e-07 +/- 2.575e-05  r1: 0.00989
        m2: 3.012e-03 +/- 7.893e-04 c2: -5.842e-06 +/- 2.700e-05  r2: -0.0404
        m:  3.038e-03 +/- 5.378e-04
        c1: 2.751e-07 +/- 2.639e-05
        c2: -5.877e-06 +/- 2.640e-05

    - run-bd24mcal01
        - fix noise run
            all but 250 or so completed
          m1: -3.352e-04 +/- 6.164e-04 c1: 1.385e-05 +/- 2.170e-05  r1: 0.00989
          m2: -1.303e-03 +/- 6.419e-04 c2: 3.397e-05 +/- 2.196e-05  r2: -0.0404
          m:  -8.046e-04 +/- 4.452e-04
          c1: 1.369e-05 +/- 2.184e-05
          c2: 3.329e-05 +/- 2.185e-05
            R: [[  4.61392474e-01  -5.63179994e-05]
                [ -2.94689880e-05   4.61878152e-01]]


    - run-bd24pcal-t01
        - trying pcal
        - this is a simp run
    - run-bd24pcal-t02
        - wider priors

- sim-bd25
    - psf size 1.5 and gal size from cosmos 23.5 sample, lognorm mean ~1.5 sigma ~1.0
      m1: 7.059e-04 +/- 7.400e-04 c1: 4.845e-05 +/- 2.564e-05
      m2: 1.854e-04 +/- 7.040e-04 c2: 5.687e-05 +/- 2.516e-05
      m:  4.375e-04 +/- 5.107e-04
      c1: 4.832e-05 +/- 2.542e-05
      c2: 5.693e-05 +/- 2.542e-05
    - so it isn't smallness of galaxies, but smallness relative to the PSF, that
        is the problem

- sim-bd26
    - low noise, no detrend
          m1: -1.433e-03 +/- 4.588e-04 c1: -1.828e-06 +/- 1.650e-05  r1: -0.13
          m2: -6.354e-04 +/- 5.440e-04 c2: 6.356e-05 +/- 1.888e-05  r2: -0.00774


- sim-bd29
    - somewhat similar to bd23
    - moffat but drawing fwhm from lognormal similar
        to DES
    - run-bd29mcal01
        200 million run
          m1: -3.933e-04 +/- 6.790e-04 c1: 5.177e-06 +/- 2.385e-05  r1: 0.0223
          m2: -3.105e-04 +/- 6.843e-04 c2: 9.410e-05 +/- 2.350e-05  r2: -0.0198
          m:  -3.528e-04 +/- 4.823e-04
          c1: 5.208e-06 +/- 2.369e-05
          c2: 9.413e-05 +/- 2.369e-05

        - max like
          m1: -4.958e-01 +/- 2.921e-04 c1: 1.145e-05 +/- 1.026e-05  r1: 0.0223
          m2: -4.955e-01 +/- 3.127e-04 c2: 1.510e-03 +/- 1.073e-05  r2: -0.0198
          m:  -4.956e-01 +/- 2.140e-04
          c1: 1.158e-05 +/- 1.051e-05
          c2: 1.510e-03 +/- 1.051e-05

- sim-bd29stars
    - run-bd29stars-mcal01
        combined with run-bd29mcal
          m1: -2.741e-04 +/- 7.027e-04 c1: 1.959e-06 +/- 2.468e-05  r1: 0.0223
          m2: -1.984e-04 +/- 7.256e-04 c2: 1.093e-04 +/- 2.491e-05  r2: -0.0198
          m:  -2.371e-04 +/- 5.052e-04
          c1: 1.988e-06 +/- 2.481e-05
          c2: 1.094e-04 +/- 2.481e-05
        maximum likelihood
          m1: -5.413e-01 +/- 2.755e-04 c1: 1.077e-05 +/- 9.674e-06  r1: 0.0223
          m2: -5.409e-01 +/- 2.992e-04 c2: 1.829e-03 +/- 1.027e-05  r2: -0.0198
          m:  -5.411e-01 +/- 2.034e-04
          c1: 1.090e-05 +/- 9.987e-06
          c2: 1.829e-03 +/- 9.987e-06

    - run-bd29stars-mcal02
        combined all run-bd29mcal01 run-bd29mcal02 run-bd29stars-mcal01 run-bd29stars-mcal02
          m1: -3.341e-04 +/- 5.013e-04 c1: 3.023e-06 +/- 1.761e-05  r1: 0.0223
          m2: -9.095e-05 +/- 5.136e-04 c2: 1.206e-04 +/- 1.763e-05  r2: -0.0198
          m:  -2.153e-04 +/- 3.590e-04
          c1: 3.116e-06 +/- 1.763e-05
          c2: 1.207e-04 +/- 1.763e-05
        max likelihood
          m1: -5.417e-01 +/- 1.975e-04 c1: 6.938e-06 +/- 6.938e-06  r1: 0.0223
          m2: -5.415e-01 +/- 2.072e-04 c2: 1.840e-03 +/- 7.114e-06  r2: -0.0198
          m:  -5.416e-01 +/- 1.432e-04
          c1: 7.004e-06 +/- 7.031e-06
          c2: 1.840e-03 +/- 7.031e-06




- run-bd29mcal02
    - Mike's rotate noise idea
          m1: -4.208e-04 +/- 6.575e-04 c1: -3.798e-06 +/- 2.309e-05  r1: 0.0223
          m2: 4.346e-05 +/- 6.863e-04 c2: 1.292e-04 +/- 2.357e-05  r2: -0.0198
          m:  -1.939e-04 +/- 4.754e-04
          c1: -3.620e-06 +/- 2.335e-05
          c2: 1.294e-04 +/- 2.335e-05
    -combined with bd29stars-mcal02
          m1: -3.942e-04 +/- 6.818e-04 c1: 4.155e-06 +/- 2.395e-05  r1: 0.0223
          m2: 1.826e-05 +/- 7.104e-04 c2: 1.321e-04 +/- 2.439e-05  r2: -0.0198
          m:  -1.926e-04 +/- 4.925e-04
          c1: 4.313e-06 +/- 2.419e-05
          c2: 1.322e-04 +/- 2.419e-05

- combined 01 and 02
      m1: -4.072e-04 +/- 4.830e-04 c1: 6.922e-07 +/- 1.696e-05  r1: 0.0223
      m2: -1.346e-04 +/- 4.892e-04 c2: 1.116e-04 +/- 1.680e-05  r2: -0.0198
      m:  -2.740e-04 +/- 3.440e-04
      c1: 7.965e-07 +/- 1.689e-05
      c2: 1.116e-04 +/- 1.689e-05

      max like
      m1: -4.958e-01 +/- 2.095e-04 c1: 4.918e-06 +/- 7.356e-06  r1: 0.0223
      m2: -4.956e-01 +/- 2.187e-04 c2: 1.521e-03 +/- 7.509e-06  r2: -0.0198
      m:  -4.957e-01 +/- 1.515e-04
      c1: 4.995e-06 +/- 7.438e-06
      c2: 1.521e-03 +/- 7.438e-06

tests of full wcs
    - run-bd29mcal-t04
        - this might be exactly reversed, so would not notice in this calculation
        - JacobianWCS(dudrow,dudcol,dvdrow,dvdcol)
              m1: -3.001e-03 +/- 6.635e-03 c1: -1.659e-04 +/- 2.330e-04  r1: 0.0223
              m2: -1.916e-02 +/- 6.880e-03 c2: -8.835e-05 +/- 2.362e-04  r2: -0.0198
              m:  -1.090e-02 +/- 4.784e-03
              c1: -1.720e-04 +/- 2.350e-04
              c2: -9.396e-05 +/- 2.349e-04
    - run-bd29mcal-t04b
        - 10 million
            m1: -3.326e-03 +/- 3.057e-03 c1: 3.015e-05 +/- 1.074e-04  r1: 0.0223
            m2: 2.848e-03 +/- 3.028e-03 c2: 1.908e-05 +/- 1.040e-04  r2: -0.0198
            m:  -3.083e-04 +/- 2.154e-03
            c1: 3.251e-05 +/- 1.058e-04
            c2: 2.123e-05 +/- 1.058e-04

    - run-bd29mcal-t05
        - Mike's suggestion
        - JacobianWCS(dudcol,dudrow,dvdcol,dvdrow)
              m1: -1.996e+00 +/- 6.683e-03 c1: -2.218e-04 +/- 2.347e-04  r1: 0.0223
              m2: 7.540e-03 +/- 6.874e-03 c2: 2.545e-04 +/- 2.360e-04  r2: -0.0198
              m:  -1.017e+00 +/- 2.292e-02
              c1: 5.452e-04 +/- 1.125e-03
              c2: 9.508e-04 +/- 1.125e-03
    -run-bd29mcal-t05c
        - after adjusting ngmix conventions.  Too noisy to be totally
            conclusive
              m1: 3.928e-03 +/- 6.799e-03 c1: 3.412e-04 +/- 2.388e-04  r1: 0.0223
              m2: -4.129e-04 +/- 6.885e-03 c2: 5.653e-05 +/- 2.364e-04  r2: -0.0198
              m:  1.807e-03 +/- 4.841e-03
              c1: 3.395e-04 +/- 2.378e-04
              c2: 5.502e-05 +/- 2.378e-04

    - run-bd29mcal-t06
        - this might be exactly reversed, so would not notice in this calculation
        - JacobianWCS(dvdcol,dvdrow,dudcol,dudrow)
              m1: -2.453e-04 +/- 6.583e-03 c1: 1.897e-05 +/- 2.312e-04  r1: 0.0223
              m2: 1.031e-04 +/- 6.735e-03 c2: 3.031e-04 +/- 2.312e-04  r2: -0.0198
              m:  -7.505e-05 +/- 4.711e-03
              c1: 1.910e-05 +/- 2.314e-04
              c2: 3.032e-04 +/- 2.314e-04
    - run-bd29mcal-t06b
        - 10 million
              m1: -2.560e-03 +/- 3.044e-03 c1: -1.382e-04 +/- 1.069e-04  r1: 0.0223
              m2: 2.633e-03 +/- 2.980e-03 c2: -9.002e-05 +/- 1.023e-04  r2: -0.0198
              m:  -2.206e-05 +/- 2.133e-03
              c1: -1.362e-04 +/- 1.047e-04
              c2: -8.822e-05 +/- 1.047e-04

    - run-bd29mcal-t07
        - testing new pixel in metacal
              m1: 3.031e-03 +/- 3.154e-03 c1: -1.737e-04 +/- 1.108e-04  r1: 0.0223
              m2: -4.077e-03 +/- 3.330e-03 c2: -6.198e-05 +/- 1.143e-04  r2: -0.0198
              m:  -4.428e-04 +/- 2.295e-03
              c1: -1.764e-04 +/- 1.127e-04
              c2: -6.445e-05 +/- 1.127e-04


- sim-bd30
    - complex wcs
- sim-bd31
    - complex wcs
    - grid of shears
    - run-bd31mcal01
        # partial run
        m1: 5.269e-04 +/- 8.985e-04 c1: -6.240e-04 +/- 3.151e-05  r1: 6.48e-18
        m2: -9.701e-04 +/- 8.712e-04 c2: 4.938e-04 +/- 3.055e-05  r2: -5.76e-18
        m:  -2.216e-04 +/- 6.264e-04
        c1: -6.240e-04 +/- 3.107e-05
        c2: 4.938e-04 +/- 3.107e-05
    - clearly something very wrong with the psf shear correction


    - run-bd31mcal-t02
        - mike's new suggestion for the pixel

- run-bd31mcal-t03
  m1: -9.362e-03 +/- 6.493e-03 c1: -6.673e-04 +/- 2.277e-04  r1: 6.48e-18
  m2: -8.912e-03 +/- 6.686e-03 c2: 6.271e-04 +/- 2.345e-04  r2: -5.76e-18
  m:  -9.137e-03 +/- 4.664e-03
  c1: -6.673e-04 +/- 2.313e-04
  c2: 6.271e-04 +/- 2.313e-04
- run-bd31mcal-t04
    - using regular fixnoise
- run-bd31mcal-t05
- run-bd31mcal-t06
    - same as t03 but now with Mike's fix in the deconvolution, where
        1/smallnumber was happening for very sheared wcs
      m1: 6.767e-03 +/- 3.122e-03 c1: -3.862e-04 +/- 1.095e-04
      m2: -2.213e-03 +/- 3.174e-03 c2: 7.835e-04 +/- 1.113e-04

 - run-bd32mcal-t01
  m1: 1.722e-03 +/- 6.500e-03 c1: -5.394e-04 +/- 2.279e-04  r1: 6.48e-18
  m2: -7.309e-04 +/- 6.849e-03 c2: 1.083e-04 +/- 2.402e-04  r2: -5.76e-18
  m:  4.954e-04 +/- 4.725e-03
  c1: -5.394e-04 +/- 2.343e-04
  c2: 1.083e-04 +/- 2.343e-04

 - run-bd33mcal-t01
  m1: 1.648e-03 +/- 6.553e-03 c1: -5.598e-04 +/- 2.298e-04  r1: 6.48e-18
  m2: 7.067e-03 +/- 6.772e-03 c2: -2.632e-04 +/- 2.375e-04  r2: -5.76e-18
  m:  4.357e-03 +/- 4.715e-03
  c1: -5.598e-04 +/- 2.339e-04
  c2: -2.632e-04 +/- 2.339e-04


- sim-bd34
    - crazy wcs like bd31
    - limiting fwhm for this crazy wcs
      m1: 6.608e-05 +/- 3.308e-03 c1: -5.381e-04 +/- 1.160e-04  r1: 6.48e-18
      m2: 3.071e-04 +/- 3.375e-03 c2: 2.330e-04 +/- 1.184e-04  r2: -5.76e-18
      m:  1.866e-04 +/- 2.365e-03
      c1: -5.381e-04 +/- 1.173e-04
      c2: 2.330e-04 +/- 1.173e-04

      still big c1

- sim-bd35
    - wcs like real DES data, with variations
    - run-bd35mcal-t01
      m1: 1.053e-03 +/- 3.003e-03 c1: -6.642e-05 +/- 1.053e-04
      m2: 2.473e-03 +/- 3.057e-03 c2: 6.086e-05 +/- 1.072e-04
      m:  1.763e-03 +/- 2.144e-03
      c1: -6.642e-05 +/- 1.063e-04
      c2: 6.086e-05 +/- 1.063e-04
    - run-bd35mcal01
      m1: -1.163e-03 +/- 6.487e-04 c1: 8.581e-06 +/- 2.275e-05  r1: 6.48e-18
      m2: 6.753e-04 +/- 6.335e-04 c2: 1.479e-04 +/- 2.222e-05  r2: -5.76e-18
      m:  -2.440e-04 +/- 4.541e-04
      c1: 8.581e-06 +/- 2.252e-05
      c2: 1.479e-04 +/- 2.252e-05
    - run-bd35mcal02
        - adding prepix psf shape, same seed
            - c1,c2 look great,but big m bias
      m1: -2.553e-03 +/- 6.680e-04 c1: 1.663e-05 +/- 2.343e-05  r1: 6.48e-18
      m2: -9.024e-04 +/- 6.456e-04 c2: 2.104e-05 +/- 2.264e-05  r2: -5.76e-18
      m:  -1.728e-03 +/- 4.652e-04
      c1: 1.663e-05 +/- 2.307e-05
      c2: 2.104e-05 +/- 2.307e-05

    - run-bd35mcal03
        - revert best I can, same seed: looks fine
      m1: -1.249e-03 +/- 6.756e-04 c1: 6.680e-06 +/- 2.369e-05  r1: 6.48e-18
      m2: 4.142e-04 +/- 6.643e-04 c2: 1.571e-04 +/- 2.330e-05  r2: -5.76e-18
      m:  -4.175e-04 +/- 4.744e-04
      c1: 6.680e-06 +/- 2.353e-05
      c2: 1.571e-04 +/- 2.353e-05
        ngmix/metacalwcs: b7b249c0197ea5af8071f44d9233a82b1fef8a10
        nsim:             b9b729dec16bbae9e143a1a8c0dbbe942c5590bf
        galsim/#720:      993044f1eb13d577e27d6f48fed23b5059ed8387

    - run-bd35mcal04
        - back to head of the branches but without prepix
        - also fixed rng

        # branch metacalwcs
        ngmix/metacalwcs fccf742687df5c59b9daa8756117a5430d04dcb4
        nsim/prepix      b9b59c2f78f5f3c9741ae386dfcd2fecf390b442
        galsim/#720      993044f1eb13d577e27d6f48fed23b5059ed8387
      m1: -2.677e-03 +/- 7.411e-04 c1: -4.449e-06 +/- 2.599e-05  r1: 6.48e-18
      m2: 2.500e-04 +/- 7.393e-04 c2: 1.405e-04 +/- 2.593e-05  r2: -5.76e-18
      m:  -1.214e-03 +/- 5.248e-04
      c1: -4.449e-06 +/- 2.603e-05
      c2: 1.405e-04 +/- 2.603e-05

- sim-bd36
    - DES wcs with -1 on dudy and dvdx, but no variation
    - run-bd36mcal-t01
        500,000
        ngmix/metacalwcs: f523769f24cba190e881569457a5d870ca198058
        nsim/prepix:      ed45fce2c18cfe0d97779dde63da0ae31b363800
        galsim/#720       993044f1eb13d577e27d6f48fed23b5059ed8387

    - run-bd36mcal-t02
        10 million
        ngmix/metacalwcs: f523769f24cba190e881569457a5d870ca198058
        nsim/prepix:      4f42b9bacdcc31c1bcbb9c277ed1b026cab11aa6
        galsim/#720       993044f1eb13d577e27d6f48fed23b5059ed8387

          m1: 2.957e-03 +/- 2.882e-03 c1: -7.198e-05 +/- 1.042e-04  r1: -0.0143
          m2: 1.200e-03 +/- 3.043e-03 c2: 2.458e-04 +/- 1.040e-04  r2: 0.0275
          m:  2.128e-03 +/- 2.094e-03
          c1: -7.155e-05 +/- 1.042e-04
          c2: 2.466e-04 +/- 1.042e-04

    - run-bd36mcal01
        - prepix
        - 100 million

          m1: 2.477e-04 +/- 8.874e-04 c1: 2.624e-05 +/- 3.209e-05  r1: -0.0143
          m2: -1.203e-03 +/- 9.065e-04 c2: 8.504e-06 +/- 3.097e-05  r2: 0.0275
          m:  -4.362e-04 +/- 6.346e-04
          c1: 2.659e-05 +/- 3.156e-05
          c2: 9.223e-06 +/- 3.157e-05

        - using two component R
          m1: 2.491e-04 +/- 8.873e-04 c1: 2.615e-05 +/- 3.208e-05  r1: -0.0143
          m2: -1.202e-03 +/- 9.064e-04 c2: 8.534e-06 +/- 3.097e-05  r2: 0.0275
          m:  -4.351e-04 +/- 6.345e-04
          c1: 2.650e-05 +/- 3.156e-05
          c2: 9.254e-06 +/- 3.156e-05

- refactored metacal, saving individual pars instead of R
    - run-bd36mcal-t03
        500,000
            ngmix/mcalrefactor: f1a8535b199be605de227fe0a4157ecbf912c89c
            nsim/prepix: fdb2ed197ebff90238344c823e24eb68af541b0a
            galsim/#720       993044f1eb13d577e27d6f48fed23b5059ed8387

            ngmix/mcalrefactor: 128d86904f068a2b85f3bbdb642685bd24829015
            nsim/prepix: d583592b53f65ab3383b4a4602e8c7ef1480d073
            galsim/#720       993044f1eb13d577e27d6f48fed23b5059ed8387
    - run-bd36mcal-t04
        10,000,000
          m1: -3.964e-03 +/- 2.889e-03 c1: -1.385e-05 +/- 1.044e-04  r1: -0.0143
          m2: -1.573e-03 +/- 3.056e-03 c2: 1.158e-04 +/- 1.044e-04  r2: 0.0275
          m:  -2.836e-03 +/- 2.101e-03
          c1: -1.443e-05 +/- 1.045e-04
          c2: 1.146e-04 +/- 1.045e-04
    - run-bd36mcal03
          100 M
            ngmix/mcalrefactor: 128d86904f068a2b85f3bbdb642685bd24829015
            nsim/prepix:        8373ee9742ba766bff1afff9de18709d8dc0a759
            galsim/#720         993044f1eb13d577e27d6f48fed23b5059ed8387

          m1: -6.561e-04 +/- 8.938e-04 c1: -3.396e-05 +/- 3.232e-05  r1: -0.0143
          m2: -3.591e-04 +/- 1.001e-03 c2: 5.506e-05 +/- 3.422e-05  r2: 0.0275
          m:  -5.160e-04 +/- 6.695e-04
          c1: -3.403e-05 +/- 3.330e-05
          c2: 5.491e-05 +/- 3.330e-05


          s/n > 7 (multiply by ~ sqrt(2) = 9.9
          m1: -1.991e-04 +/- 8.999e-04 c1: -4.301e-05 +/- 3.254e-05  r1: -0.0143
          m2: 3.008e-04 +/- 1.028e-03 c2: 3.388e-05 +/- 3.511e-05  r2: 0.0275
          m:  3.667e-05 +/- 6.810e-04
          c1: -4.314e-05 +/- 3.387e-05
          c2: 3.363e-05 +/- 3.387e-05

          s/n > 10 (*sqrt(2) = 14)
          m1: 1.401e-05 +/- 9.466e-04 c1: -1.787e-05 +/- 3.423e-05  r1: -0.0143
          m2: 1.447e-04 +/- 1.065e-03 c2: 4.009e-05 +/- 3.638e-05  r2: 0.0275
          m:  7.562e-05 +/- 7.105e-04
          c1: -1.790e-05 +/- 3.534e-05
          c2: 4.003e-05 +/- 3.534e-05
          without correction
              m1: 3.835e-03 +/- 9.503e-04 c1: -1.793e-05 +/- 3.436e-05  r1: -0.0143
              m2: 3.728e-03 +/- 1.069e-03 c2: 1.079e-04 +/- 3.651e-05  r2: 0.0275
              m:  3.784e-03 +/- 7.132e-04
              c1: -1.791e-05 +/- 3.547e-05
              c2: 1.080e-04 +/- 3.548e-05

    - run-bd36mcal04
          300 M
            ngmix/mcalrefactor: 128d86904f068a2b85f3bbdb642685bd24829015
            nsim/prepix:        5f274ad1f22430334bfbd9a9c3fc4769f008f2b3
            galsim/#720         993044f1eb13d577e27d6f48fed23b5059ed8387

    - run-bd36mcal05
          400 M
            ngmix/master:   8fa011902b3d60c286e329a4ad1836528ab9a3c4
            nsim/master:    6e4c7e5547afc42a164dc2a2c0ab5d30a0b5f9fd
            galsim/#720     993044f1eb13d577e27d6f48fed23b5059ed8387

         - combined 03/04/05
              m1: -7.575e-04 +/- 3.208e-04 c1: 4.936e-06 +/- 1.160e-05
              m2: -1.020e-03 +/- 3.394e-04 c2: 2.506e-05 +/- 1.160e-05
              m:  -8.813e-04 +/- 2.333e-04
              c1: 5.000e-06 +/- 1.161e-05
              c2: 2.519e-05 +/- 1.161e-05 

        - s/n > 3.5 (4.9)
            without correction
              m1: -6.023e-04 +/- 3.248e-04 c1: 2.578e-06 +/- 1.174e-05  r1: -0.0143
              m2: -8.845e-04 +/- 3.438e-04 c2: 2.961e-05 +/- 1.175e-05  r2: 0.0275
              m:  -7.354e-04 +/- 2.363e-04
              c1: 2.647e-06 +/- 1.175e-05
              c2: 2.975e-05 +/- 1.175e-05
            with correction
              m1: -6.650e-04 +/- 3.248e-04 c1: 2.578e-06 +/- 1.174e-05  r1: -0.0143
              m2: -9.475e-04 +/- 3.437e-04 c2: 2.793e-05 +/- 1.174e-05  r2: 0.0275
              m:  -7.982e-04 +/- 2.363e-04
              c1: 2.647e-06 +/- 1.175e-05
              c2: 2.807e-05 +/- 1.175e-05
        - s/n > 5.0 (7.1)
            without correction
              m1: 2.917e-04 +/- 3.250e-04 c1: 4.584e-06 +/- 1.175e-05  r1: -0.0143
              m2: 5.956e-05 +/- 3.465e-04 c2: 4.811e-05 +/- 1.184e-05  r2: 0.0275
              m:  1.822e-04 +/- 2.373e-04
              c1: 4.641e-06 +/- 1.180e-05
              c2: 4.823e-05 +/- 1.180e-05
            with correction
              m1: -4.544e-04 +/- 3.247e-04 c1: 4.578e-06 +/- 1.174e-05  r1: -0.0143
              m2: -6.748e-04 +/- 3.462e-04 c2: 2.554e-05 +/- 1.183e-05  r2: 0.0275
              m:  -5.583e-04 +/- 2.371e-04
              c1: 4.632e-06 +/- 1.179e-05
              c2: 2.565e-05 +/- 1.179e-05
        - s/n > 7 (9.9)
            without correction
              m1: 2.742e-03 +/- 3.265e-04 c1: -6.660e-06 +/- 1.181e-05  r1: -0.0143
              m2: 2.700e-03 +/- 3.482e-04 c2: 1.059e-04 +/- 1.190e-05  r2: 0.0275
              m:  2.722e-03 +/- 2.384e-04
              c1: -6.650e-06 +/- 1.186e-05
              c2: 1.059e-04 +/- 1.186e-05
            with correction
              m1: 1.014e-05 +/- 3.256e-04 c1: -6.649e-06 +/- 1.177e-05  r1: -0.0143
              m2: -1.312e-04 +/- 3.472e-04 c2: 2.845e-05 +/- 1.186e-05  r2: 0.0275
              m:  -5.650e-05 +/- 2.377e-04
              c1: -6.615e-06 +/- 1.182e-05
              c2: 2.852e-05 +/- 1.183e-05
        - s/n > 10 (14.1)
            without correction
              m1: 3.443e-03 +/- 3.411e-04 c1: 5.874e-08 +/- 1.233e-05  r1: -0.0143
              m2: 3.725e-03 +/- 3.610e-04 c2: 9.555e-05 +/- 1.233e-05  r2: 0.0275
              m:  3.576e-03 +/- 2.481e-04
              c1: -9.751e-09 +/- 1.234e-05
              c2: 9.541e-05 +/- 1.234e-05
            with correction
              m1: -3.195e-04 +/- 3.398e-04 c1: 5.733e-08 +/- 1.229e-05  r1: -0.0143
              m2: -2.374e-05 +/- 3.596e-04 c2: 2.778e-05 +/- 1.229e-05  r2: 0.0275
              m:  -1.800e-04 +/- 2.472e-04
              c1: -1.511e-08 +/- 1.229e-05
              c2: 2.763e-05 +/- 1.230e-05
        - s/n > 14 (19.8)
            without correction
              m1: 2.158e-03 +/- 3.608e-04 c1: 4.646e-06 +/- 1.304e-05  r1: -0.0143
              m2: 2.714e-03 +/- 3.768e-04 c2: 6.479e-05 +/- 1.288e-05  r2: 0.0275
              m:  2.420e-03 +/- 2.608e-04
              c1: 4.510e-06 +/- 1.297e-05
              c2: 6.451e-05 +/- 1.297e-05
            with correction
              m1: -4.041e-04 +/- 3.599e-04 c1: 4.633e-06 +/- 1.301e-05  r1: -0.0143
              m2: 1.049e-04 +/- 3.759e-04 c2: 2.322e-05 +/- 1.284e-05  r2: 0.0275
              m:  -1.641e-04 +/- 2.601e-04
              c1: 4.508e-06 +/- 1.294e-05
              c2: 2.296e-05 +/- 1.294e-05

        -run-bd36mcal06
            - dilate additional 1/dim
            100M

            ngmix/dilate:   cf920178222249cbbc0c4bfdcdb8e466dc76ed6e
            nsim/master:    f9f24fef016841bea50352fb58245bbfb3a28cbb
            galsim/#720     993044f1eb13d577e27d6f48fed23b5059ed8387
        -run-bd36mcal07
            ngmix/dilate:   cf920178222249cbbc0c4bfdcdb8e466dc76ed6e
            nsim/master:    5d35ead4b8df6b6adb87ea821db7fd11dac24635
            galsim/#720     993044f1eb13d577e27d6f48fed23b5059ed8387


         combined 06/07 so far
          m1: -1.141e-03 +/- 6.080e-04 c1: 1.697e-05 +/- 2.198e-05  r1: -0.0143
          m2: -9.139e-04 +/- 6.561e-04 c2: 2.542e-05 +/- 2.242e-05  r2: 0.0275
          m:  -1.034e-03 +/- 4.466e-04
          c1: 1.691e-05 +/- 2.221e-05
          c2: 2.531e-05 +/- 2.222e-05

          s/n > 7 (9.9)
            without correction
              m1: 1.515e-03 +/- 6.260e-04 c1: 1.980e-05 +/- 2.263e-05  r1: -0.0143
              m2: 2.694e-03 +/- 6.742e-04 c2: 1.052e-04 +/- 2.304e-05  r2: 0.0275
              m:  2.071e-03 +/- 4.596e-04
              c1: 1.951e-05 +/- 2.286e-05
              c2: 1.046e-04 +/- 2.286e-05

            with correction
              m1: -9.070e-04 +/- 6.245e-04 c1: 1.973e-05 +/- 2.258e-05  r1: -0.0143
              m2: 3.425e-04 +/- 6.726e-04 c2: 2.861e-05 +/- 2.298e-05  r2: 0.0275
              m:  -3.177e-04 +/- 4.585e-04
              c1: 1.943e-05 +/- 2.280e-05
              c2: 2.799e-05 +/- 2.281e-05


    - run-bd36mcal-t05
        - keeping prepix psf and that shape
            ngmix/keep_prepix: 40ba9d5e6afce5422e47c26aaf209cd21e74a5fc
            nsim/master: c1f42392907bfbf2c56e0bf953f9c7ea0f909325
            galsim/#720: e11369f15e70d8b8b5e3ef8024bbc0f5f1dd0113
    - run-bd36mcal08
        - keeping prepix psf and that shape
            ngmix/keep_prepix: 40ba9d5e6afce5422e47c26aaf209cd21e74a5fc
            nsim/master: 252230e671ab9a33d271bb481114542a0e2eecfd
            galsim/#720: 993044f1eb13d577e27d6f48fed23b5059ed8387
    - run-bd36mcal09
            ngmix/keep_prepix: 40ba9d5e6afce5422e47c26aaf209cd21e74a5fc
            nsim/master: 252230e671ab9a33d271bb481114542a0e2eecfd
            galsim/#720: 993044f1eb13d577e27d6f48fed23b5059ed8387

    - combined 08/09

        R: [ 0.41263295  0.41308002]
        Rpsf: [ 0.18915615  0.18881478]
          m1: -7.985e-04 +/- 4.461e-04 c1: 4.194e-06 +/- 1.613e-05  r1: -0.0143
          m2: -1.180e-03 +/- 4.879e-04 c2: 1.130e-06 +/- 1.667e-05  r2: 0.0275
          m:  -9.783e-04 +/- 3.300e-04
          c1: 4.287e-06 +/- 1.641e-05
          c2: 1.319e-06 +/- 1.641e-05

        selections

        s/n > 5 (7.1)
            R: [ 0.41997298  0.4204284 ]
            Rpsf: [ 0.19386113  0.19350501]

            Rsel: [ 0.00026776  0.00025949]
            Rpsf_sel: [ 0.00134689  0.00133816]
            without correction
              m1: 1.251e-04 +/- 4.481e-04 c1: 5.300e-06 +/- 1.620e-05  r1: -0.0143
              m2: -2.727e-04 +/- 4.899e-04 c2: 1.989e-05 +/- 1.674e-05  r2: 0.0275
              m:  -6.253e-05 +/- 3.314e-04
              c1: 5.397e-06 +/- 1.648e-05
              c2: 2.008e-05 +/- 1.649e-05
            with correction
              m1: -5.123e-04 +/- 4.478e-04 c1: 5.302e-06 +/- 1.619e-05  r1: -0.0143
              m2: -8.894e-04 +/- 4.896e-04 c2: -2.424e-06 +/- 1.673e-05  r2: 0.0275
              m:  -6.901e-04 +/- 3.312e-04
              c1: 5.394e-06 +/- 1.647e-05
              c2: -2.237e-06 +/- 1.648e-05
        s/n > 7 (9.9)
            R: [ 0.45748975  0.45799785]
            Rpsf: [ 0.21168433  0.21127212]

            Rsel: [ 0.00112809  0.00113626]
            Rpsf_sel: [ 0.00482381  0.0048189 ]
            without correction
              m1: 2.510e-03 +/- 4.518e-04 c1: -5.931e-06 +/- 1.633e-05  r1: -0.0143
              m2: 1.940e-03 +/- 5.001e-04 c2: 7.178e-05 +/- 1.709e-05  r2: 0.0275
              m:  2.241e-03 +/- 3.363e-04
              c1: -5.792e-06 +/- 1.673e-05
              c2: 7.206e-05 +/- 1.673e-05
            with correction
              m1: 4.330e-05 +/- 4.506e-04 c1: -5.899e-06 +/- 1.629e-05  r1: -0.0143
              m2: -5.397e-04 +/- 4.989e-04 c2: -1.974e-06 +/- 1.704e-05  r2: 0.0275
              m:  -2.317e-04 +/- 3.355e-04
              c1: -5.756e-06 +/- 1.669e-05
              c2: -1.685e-06 +/- 1.669e-05
        s/n > 10 (14.1)
            R: [ 0.53680552  0.53739111]
            Rpsf: [ 0.23321205  0.23269319]

            Rsel: [ 0.00170984  0.00175894]
            Rpsf_sel: [ 0.00480542  0.00479401]
            without correction
              m1: 3.612e-03 +/- 4.696e-04 c1: -5.731e-06 +/- 1.698e-05  r1: -0.0143
              m2: 3.072e-03 +/- 5.200e-04 c2: 5.808e-05 +/- 1.777e-05  r2: 0.0275
              m:  3.357e-03 +/- 3.496e-04
              c1: -5.599e-06 +/- 1.739e-05
              c2: 5.835e-05 +/- 1.739e-05
            with correction
              m1: 4.254e-04 +/- 4.681e-04 c1: -5.705e-06 +/- 1.693e-05  r1: -0.0143
              m2: -2.009e-04 +/- 5.183e-04 c2: -4.443e-06 +/- 1.771e-05  r2: 0.0275
              m:  1.301e-04 +/- 3.485e-04
              c1: -5.552e-06 +/- 1.734e-05
              c2: -4.132e-06 +/- 1.734e-05

        s/n > 14 (19.8)
            R: [ 0.60200558  0.60269422]
            Rpsf: [ 0.24430521  0.24369149]

            Rsel: [ 0.00126953  0.00131957]
            Rpsf_sel: [ 0.00341135  0.00338747]
            without correction
              m1: 2.312e-03 +/- 4.817e-04 c1: 7.853e-07 +/- 1.742e-05  r1: -0.0143
              m2: 2.548e-03 +/- 5.246e-04 c2: 4.234e-05 +/- 1.792e-05  r2: 0.0275
              m:  2.423e-03 +/- 3.555e-04
              c1: 7.278e-07 +/- 1.768e-05
              c2: 4.222e-05 +/- 1.769e-05
            with correction
              m1: 2.028e-04 +/- 4.807e-04 c1: 7.928e-07 +/- 1.738e-05  r1: -0.0143
              m2: 3.577e-04 +/- 5.234e-04 c2: 2.938e-06 +/- 1.788e-05  r2: 0.0275
              m:  2.758e-04 +/- 3.548e-04
              c1: 7.549e-07 +/- 1.764e-05
              c2: 2.861e-06 +/- 1.765e-05

    - combined 3,4,5,8,9 1200M

        R: [ 0.41585284  0.41626658]
        Rpsf: [ 0.1953022  0.1950376]
          m1: -7.705e-04 +/- 2.533e-04 c1: 4.699e-06 +/- 9.158e-06  r1: -0.0143
          m2: -1.073e-03 +/- 2.806e-04 c2: 1.715e-05 +/- 9.588e-06  r2: 0.0275
          m:  -9.132e-04 +/- 1.886e-04
          c1: 4.773e-06 +/- 9.382e-06
          c2: 1.730e-05 +/- 9.383e-06

        s/n > 5 (7.1)   0.97
            R: [ 0.42300142  0.42341702]
            Rpsf: [ 0.19992092  0.1996469 ]
            Rsel: [ 0.00030157  0.00029417]
            Rpsf_sel: [ 0.00135974  0.00135685]
            without correction
              m1: 2.016e-04 +/- 2.541e-04 c1: 6.112e-06 +/- 9.186e-06  r1: -0.0143
              m2: -8.053e-05 +/- 2.819e-04 c2: 3.726e-05 +/- 9.633e-06  r2: 0.0275
              m:  6.853e-05 +/- 1.894e-04
              c1: 6.181e-06 +/- 9.419e-06
              c2: 3.740e-05 +/- 9.420e-06
            with correction
              m1: -5.110e-04 +/- 2.539e-04 c1: 6.108e-06 +/- 9.179e-06  r1: -0.0143
              m2: -7.747e-04 +/- 2.818e-04 c2: 1.478e-05 +/- 9.627e-06  r2: 0.0275
              m:  -6.354e-04 +/- 1.892e-04
              c1: 6.172e-06 +/- 9.412e-06
              c2: 1.491e-05 +/- 9.413e-06

        s/n > 7 (9.9)   0.83
            R: [ 0.45944856  0.45991696]
            Rpsf: [ 0.2173346  0.2170133]
            Rsel: [ 0.00121558  0.00124866]
            Rpsf_sel: [ 0.00502684  0.00499741]
            without correction
              m1: 2.616e-03 +/- 2.587e-04 c1: -5.312e-06 +/- 9.354e-06  r1: -0.0143
              m2: 2.424e-03 +/- 2.883e-04 c2: 9.344e-05 +/- 9.850e-06  r2: 0.0275
              m:  2.526e-03 +/- 1.932e-04
              c1: -5.265e-06 +/- 9.611e-06
              c2: 9.354e-05 +/- 9.613e-06
            with correction
              m1: -2.945e-05 +/- 2.580e-04 c1: -5.297e-06 +/- 9.329e-06  r1: -0.0143
              m2: -2.905e-04 +/- 2.875e-04 c2: 1.722e-05 +/- 9.823e-06  r2: 0.0275
              m:  -1.526e-04 +/- 1.927e-04
              c1: -5.233e-06 +/- 9.586e-06
              c2: 1.735e-05 +/- 9.587e-06
        s/n > 10 (14.1)   0.60
            R: [ 0.53788513  0.53842082]
            Rpsf: [ 0.23816384  0.2377501 ]
            Rsel: [ 0.00191993  0.0019345 ]
            Rpsf_sel: [ 0.005073    0.00506797]
            without correction
              m1: 3.447e-03 +/- 2.674e-04 c1: -3.316e-07 +/- 9.668e-06  r1: -0.0143
              m2: 3.513e-03 +/- 2.983e-04 c2: 8.152e-05 +/- 1.019e-05  r2: 0.0275
              m:  3.478e-03 +/- 1.998e-04
              c1: -3.475e-07 +/- 9.939e-06
              c2: 8.149e-05 +/- 9.940e-06

            with correction
              m1: -1.223e-04 +/- 2.664e-04 c1: -3.274e-07 +/- 9.634e-06  r1: -0.0143
              m2: -8.042e-05 +/- 2.972e-04 c2: 1.549e-05 +/- 1.016e-05  r2: 0.0275
              m:  -1.025e-04 +/- 1.991e-04
              c1: -3.376e-07 +/- 9.904e-06
              c2: 1.547e-05 +/- 9.905e-06

        s/n > 13 (18.4)  0.46
            R: [ 0.58981424  0.5904208 ]
            Rpsf: [ 0.2469236  0.246439 ]
            Rsel: [ 0.00157791  0.00161063]
            Rpsf_sel: [ 0.00380999  0.00379327]
            without correction
              m1: 2.374e-03 +/- 2.750e-04 c1: 4.310e-06 +/- 9.941e-06  r1: -0.0143
              m2: 2.873e-03 +/- 3.048e-04 c2: 5.819e-05 +/- 1.042e-05  r2: 0.0275
              m:  2.609e-03 +/- 2.049e-04
              c1: 4.189e-06 +/- 1.019e-05
              c2: 5.794e-05 +/- 1.019e-05
            with correction
              m1: -3.004e-04 +/- 2.742e-04 c1: 4.301e-06 +/- 9.915e-06  r1: -0.0143
              m2: 1.445e-04 +/- 3.040e-04 c2: 1.312e-05 +/- 1.039e-05  r2: 0.0275
              m:  -9.059e-05 +/- 2.043e-04
              c1: 4.192e-06 +/- 1.016e-05
              c2: 1.290e-05 +/- 1.016e-05

        s/n > 16 (22.6) 0.39
            R: [ 0.6247135   0.62538052]
            Rpsf: [ 0.25055988  0.25001014]
            Rsel: [ 0.00117735  0.001183  ]
            Rpsf_sel: [ 0.00316194  0.00313398]
            without correction
              m1: 1.571e-03 +/- 2.848e-04 c1: 7.975e-06 +/- 1.030e-05  r1: -0.0143
              m2: 2.158e-03 +/- 3.223e-04 c2: 4.940e-05 +/- 1.101e-05  r2: 0.0275
              m:  1.848e-03 +/- 2.146e-04
              c1: 7.832e-06 +/- 1.067e-05
              c2: 4.911e-05 +/- 1.067e-05
            with correction
              m1: -3.131e-04 +/- 2.843e-04 c1: 7.964e-06 +/- 1.028e-05  r1: -0.0143
              m2: 2.658e-04 +/- 3.217e-04 c2: 1.425e-05 +/- 1.099e-05  r2: 0.0275
              m:  -4.009e-05 +/- 2.142e-04
              c1: 7.822e-06 +/- 1.065e-05
              c2: 1.396e-05 +/- 1.065e-05

        s/n > 19 (26.9)  0.33
            R: [ 0.64969728  0.65038296]
            Rpsf: [ 0.25178775  0.25119776]
            Rsel: [ 0.00077364  0.00079763]
            Rpsf_sel: [ 0.00276496  0.00274198]
            without correction
              m1: 6.086e-04 +/- 2.953e-04 c1: 3.426e-06 +/- 1.068e-05  r1: -0.0143
              m2: 1.388e-03 +/- 3.303e-04 c2: 4.224e-05 +/- 1.129e-05  r2: 0.0275
              m:  9.761e-04 +/- 2.212e-04
              c1: 3.236e-06 +/- 1.100e-05
              c2: 4.185e-05 +/- 1.100e-05
            with correction
              m1: -5.816e-04 +/- 2.949e-04 c1: 3.425e-06 +/- 1.066e-05  r1: -0.0143
              m2: 1.612e-04 +/- 3.299e-04 c2: 1.268e-05 +/- 1.127e-05  r2: 0.0275
              m:  -2.313e-04 +/- 2.209e-04
              c1: 3.243e-06 +/- 1.099e-05
              c2: 1.231e-05 +/- 1.099e-05


          # 4 ranges

            7 < snr < 16  0.44
            ------------------

            mean: 10.3982

            R: [ 0.31493787  0.31523259]
            Rpsf: [ 0.18828181  0.18816025]
            Rsel: [ 0.00124903  0.00130613]
            Rpsf_sel: [ 0.00665751  0.00662684]
            without correction
              m1: 4.444e-03 +/- 5.098e-04 c1: -2.948e-05 +/- 1.843e-05  r1: -0.0143
              m2: 2.816e-03 +/- 5.694e-04 c2: 1.694e-04 +/- 1.946e-05  r2: 0.0275
              m:  3.676e-03 +/- 3.817e-04
              c1: -2.908e-05 +/- 1.898e-05
              c2: 1.702e-04 +/- 1.899e-05
            with correction
              m1: 4.759e-04 +/- 5.078e-04 c1: -2.937e-05 +/- 1.836e-05  r1: -0.0143
              m2: -1.323e-03 +/- 5.671e-04 c2: 2.190e-05 +/- 1.938e-05  r2: 0.0275
              m:  -3.725e-04 +/- 3.802e-04
              c1: -2.893e-05 +/- 1.891e-05
              c2: 2.279e-05 +/- 1.891e-05

            16 < snr < 37  0.22
            ------------------

            mean: 24.0881

            R: [ 0.55292432  0.55347453]
            Rpsf: [ 0.25176904  0.25137178]
            Rsel: [ 0.00274647  0.00273931]
            Rpsf_sel: [ 0.00389896  0.0038854 ]
            without correction
              m1: 4.637e-03 +/- 4.474e-04 c1: 1.219e-05 +/- 1.618e-05  r1: -0.0143
              m2: 4.453e-03 +/- 4.984e-04 c2: 7.033e-05 +/- 1.703e-05  r2: 0.0275
              m:  4.550e-03 +/- 3.341e-04
              c1: 1.223e-05 +/- 1.662e-05
              c2: 7.042e-05 +/- 1.662e-05
            with correction
              m1: -3.288e-04 +/- 4.452e-04 c1: 1.214e-05 +/- 1.610e-05  r1: -0.0143
              m2: -4.940e-04 +/- 4.960e-04 c2: 2.102e-05 +/- 1.695e-05  r2: 0.0275
              m:  -4.067e-04 +/- 3.325e-04
              c1: 1.218e-05 +/- 1.654e-05
              c2: 2.110e-05 +/- 1.654e-05


            37 < snr < 87 0.10
            ------------------

            mean: 55.41

            R: [ 0.69248894  0.69323276]
            Rpsf: [ 0.25283065  0.25216233]
            Rsel: [  1.99787337e-05   7.51765451e-05]
            Rpsf_sel: [ 0.00214137  0.00207139]
            without correction
              m1: 1.681e-05 +/- 5.112e-04 c1: -1.261e-05 +/- 1.848e-05  r1: -0.0143
              m2: 4.253e-04 +/- 5.406e-04 c2: 4.301e-06 +/- 1.847e-05  r2: 0.0275
              m:  2.095e-04 +/- 3.717e-04
              c1: -1.271e-05 +/- 1.849e-05
              c2: 4.098e-06 +/- 1.849e-05
            with correction
              m1: -1.209e-05 +/- 5.112e-04 c1: -1.261e-05 +/- 1.848e-05  r1: -0.0143
              m2: 3.166e-04 +/- 5.405e-04 c2: -1.664e-05 +/- 1.847e-05  r2: 0.0275
              m:  1.429e-04 +/- 3.717e-04
              c1: -1.269e-05 +/- 1.849e-05
              c2: -1.680e-05 +/- 1.849e-05

            87 < snr < 200 0.04
            ------------------

            mean: 127.428

            R: [ 0.74279904  0.74371095]
            Rpsf: [ 0.24467445  0.24384656]
            Rsel: [-0.00161131 -0.00169201]
            Rpsf_sel: [ 0.00217835  0.00217786]
            without correction
              m1: -2.913e-03 +/- 7.044e-04 c1: -3.709e-06 +/- 2.547e-05  r1: -0.0143
              m2: -9.734e-04 +/- 7.633e-04 c2: 8.229e-05 +/- 2.608e-05  r2: 0.0275
              m:  -1.998e-03 +/- 5.190e-04
              c1: -4.183e-06 +/- 2.581e-05
              c2: 8.133e-05 +/- 2.582e-05
            with correction
              m1: -7.455e-04 +/- 7.060e-04 c1: -3.712e-06 +/- 2.553e-05  r1: -0.0143
              m2: 1.305e-03 +/- 7.650e-04 c2: 6.192e-05 +/- 2.614e-05  r2: 0.0275
              m:  2.215e-04 +/- 5.202e-04
              c1: -4.212e-06 +/- 2.587e-05
              c2: 6.090e-05 +/- 2.588e-05


          # 3 ranges
          7 < s/n < 10.84 0.27
            R: [ 0.26768178  0.26797514]
            Rpsf: [ 0.16906917  0.16897402]
            Rsel: [ -7.15111132e-05  -7.19682243e-05]
            Rpsf_sel: [ 0.00582695  0.00579787]
            without correction
              m1: 3.471e-04 +/- 7.291e-04 c1: -3.848e-05 +/- 2.636e-05  r1: -0.0143
              m2: -1.179e-03 +/- 8.267e-04 c2: 1.776e-04 +/- 2.825e-05  r2: 0.0275
              m:  -3.727e-04 +/- 5.498e-04
              c1: -3.811e-05 +/- 2.735e-05
              c2: 1.784e-04 +/- 2.735e-05
            with correction
              m1: 6.144e-04 +/- 7.292e-04 c1: -3.850e-05 +/- 2.637e-05  r1: -0.0143
              m2: -9.122e-04 +/- 8.269e-04 c2: 2.590e-05 +/- 2.825e-05  r2: 0.0275
              m:  -1.055e-04 +/- 5.500e-04
              c1: -3.813e-05 +/- 2.735e-05
              c2: 2.666e-05 +/- 2.736e-05

          10.84 < s/n < 22.73  0.28
            R: [ 0.43752933  0.43791064]
            Rpsf: [ 0.23084     0.23060314]
            Rsel: [ 0.00338639  0.00348928]
            Rpsf_sel: [ 0.00677821  0.00671319]
            without correction
              m1: 7.966e-03 +/- 4.983e-04 c1: -3.096e-06 +/- 1.802e-05  r1: -0.0143
              m2: 6.849e-03 +/- 5.482e-04 c2: 1.155e-04 +/- 1.873e-05  r2: 0.0275
              m:  7.439e-03 +/- 3.699e-04
              c1: -2.824e-06 +/- 1.840e-05
              c2: 1.161e-04 +/- 1.840e-05
            with correction
              m1: 2.235e-04 +/- 4.945e-04 c1: -3.065e-06 +/- 1.788e-05  r1: -0.0143
              m2: -1.111e-03 +/- 5.439e-04 c2: 8.000e-06 +/- 1.858e-05  r2: 0.0275
              m:  -4.058e-04 +/- 3.671e-04
              c1: -2.740e-06 +/- 1.826e-05
              c2: 8.662e-06 +/- 1.826e-05

          s/n > 22.73
            R: [ 0.67189152  0.67262057]
            Rpsf: [ 0.25187801  0.25124789]
            Rsel: [ 0.00033528  0.00033212]
            Rpsf_sel: [ 0.00248891  0.00249481]
            without correction
              m1: 5.719e-05 +/- 3.171e-04 c1: 5.653e-06 +/- 1.147e-05  r1: -0.0143
              m2: 9.372e-04 +/- 3.416e-04 c2: 4.604e-05 +/- 1.167e-05  r2: 0.0275
              m:  4.722e-04 +/- 2.329e-04
              c1: 5.438e-06 +/- 1.159e-05
              c2: 4.560e-05 +/- 1.159e-05
            with correction
              m1: -4.417e-04 +/- 3.170e-04 c1: 5.652e-06 +/- 1.146e-05  r1: -0.0143
              m2: 4.432e-04 +/- 3.414e-04 c2: 2.003e-05 +/- 1.166e-05  r2: 0.0275
              m:  -2.441e-05 +/- 2.328e-04
              c1: 5.436e-06 +/- 1.158e-05
              c2: 1.960e-05 +/- 1.158e-05

        - max like
          m1: -4.958e-01 +/- 1.107e-04 c1: -2.130e-06 +/- 4.004e-06  r1: -0.0143
          m2: -4.956e-01 +/- 1.227e-04 c2: 1.516e-03 +/- 4.191e-06  r2: 0.0275
          m:  -4.957e-01 +/- 8.248e-05
          c1: -2.175e-06 +/- 4.102e-06
          c2: 1.516e-03 +/- 4.103e-06

    - run-bd36stars-mcal01
        nsim/master:       5d29414fbd835db70fbd14cca12aa17eba925431
        ngmix/keep_prepix: 8537b9dacde34d75907cb3dc0c92022384470a9f
        galsim/#720:       993044f1eb13d577e27d6f48fed23b5059ed8387
    - run-bd36stars-mcal02
        nsim/master:       5d29414fbd835db70fbd14cca12aa17eba925431
        ngmix/keep_prepix: 8537b9dacde34d75907cb3dc0c92022384470a9f
        galsim/#720:       993044f1eb13d577e27d6f48fed23b5059ed8387

    - run-bd36stars-mcal03
        nsim/master:       aff87d126ee0f019a64a55e2e305ae9969390127
        ngmix/keep_prepix: 60ce3ba996ea220a5c781e078611262883faf7e5
        galsim/#720:       64192dfbf2bd7f9a616cf6dae7cf7476adace9c4

        - combined with above and s/n > 7
        without correction
          m1: 2.599e-03 +/- 2.694e-04 c1: -8.826e-06 +/- 9.740e-06  r1: -0.0143
          m2: 2.464e-03 +/- 3.010e-04 c2: 9.445e-05 +/- 1.028e-05  r2: 0.0275
          m:  2.535e-03 +/- 2.015e-04
          c1: -8.793e-06 +/- 1.002e-05
          c2: 9.451e-05 +/- 1.002e-05

        with correction
          m1: -4.819e-05 +/- 2.687e-04 c1: -8.801e-06 +/- 9.714e-06  r1: -0.0143
          m2: -2.483e-04 +/- 3.001e-04 c2: 1.439e-05 +/- 1.025e-05  r2: 0.0275
          m:  -1.425e-04 +/- 2.009e-04
          c1: -8.753e-06 +/- 9.994e-06
          c2: 1.449e-05 +/- 9.996e-06

        max like

          m1: -5.417e-01 +/- 1.043e-04 c1: -3.359e-06 +/- 3.772e-06  r1: -0.0143
          m2: -5.415e-01 +/- 1.154e-04 c2: 1.833e-03 +/- 3.944e-06  r2: 0.0275
          m:  -5.416e-01 +/- 7.766e-05
          c1: -3.411e-06 +/- 3.863e-06
          c2: 1.833e-03 +/- 3.863e-06

- sim-bdmask05
    - goes with sim-bd36
- run-bdmask05mcal01
    - 200 M
    038830db2c03db9601bf0c9cd53afaee2f3dcb0d
        nsim/master:       038830db2c03db9601bf0c9cd53afaee2f3dcb0d
        ngmix/keep_prepix: 60ce3ba996ea220a5c781e078611262883faf7e5
        galsim/#720:       64192dfbf2bd7f9a616cf6dae7cf7476adace9c4
    s/n > 7
          m1: 2.651e-03 +/- 6.507e-04 c1: 3.010e-05 +/- 2.282e-05  r1: 6.48e-18
          m2: 2.371e-03 +/- 6.722e-04 c2: 1.309e-04 +/- 2.357e-05  r2: -5.76e-18
          m:  2.511e-03 +/- 4.681e-04
          c1: 3.010e-05 +/- 2.322e-05
          c2: 1.309e-04 +/- 2.322e-05

        with correction
          m1: 2.406e-04 +/- 6.491e-04 c1: 3.001e-05 +/- 2.276e-05  r1: 6.48e-18
          m2: -1.113e-05 +/- 6.706e-04 c2: 5.613e-05 +/- 2.352e-05  r2: -5.76e-18
          m:  1.147e-04 +/- 4.670e-04
          c1: 3.001e-05 +/- 2.316e-05
          c2: 5.613e-05 +/- 2.316e-05
    no cut
          m1: -9.776e-04 +/- 6.311e-04 c1: 3.052e-05 +/- 2.213e-05  r1: 6.48e-18
          m2: -7.083e-04 +/- 6.675e-04 c2: 7.304e-05 +/- 2.341e-05  r2: -5.76e-18
          m:  -8.430e-04 +/- 4.597e-04
          c1: 3.052e-05 +/- 2.280e-05
          c2: 7.304e-05 +/- 2.280e-05

- run-bd36mcal-rnd01
    - using round gaussian psf for final fits, dilated by factor
        of 1.2

    - 174/200 done
        note R is 0.36 compared to 0.42
        R: [ 0.35943831  0.359544  ]
        Rpsf: [ 0.  0.]
          m1: -1.410e-03 +/- 6.569e-04 c1: 2.990e-05 +/- 2.375e-05  r1: -0.0143
          m2: -7.128e-04 +/- 7.163e-04 c2: -2.595e-06 +/- 2.447e-05  r2: 0.0275
          m:  -1.081e-03 +/- 4.852e-04
          c1: 2.973e-05 +/- 2.413e-05
          c2: -2.941e-06 +/- 2.414e-05

        - s/n > 10
        R: [ 0.53439255  0.53451936]
        Rpsf: [ 0.  0.]

        Rsel: [-0.00024671 -0.00024456]
        Rpsf_sel: [ 0.  0.]
        without correction
          m1: -1.078e-03 +/- 7.272e-04 c1: 5.488e-05 +/- 2.629e-05  r1: -0.0143
          m2: -7.512e-05 +/- 7.484e-04 c2: 1.825e-05 +/- 2.557e-05  r2: 0.0275
          m:  -6.053e-04 +/- 5.218e-04
          c1: 5.464e-05 +/- 2.595e-05
          c2: 1.775e-05 +/- 2.596e-05

        with correction
          m1: -6.170e-04 +/- 7.275e-04 c1: 5.491e-05 +/- 2.630e-05  r1: -0.0143
          m2: 3.826e-04 +/- 7.488e-04 c2: 1.825e-05 +/- 2.558e-05  r2: 0.0275
          m:  -1.456e-04 +/- 5.221e-04
          c1: 5.466e-05 +/- 2.597e-05
          c2: 1.776e-05 +/- 2.597e-05



- sim-bd37
    - actual r50 drawn from cosmos sample used in great3
    - run-bd37zmax-lownoise01
        nsim/master: e863a1d2f2ee817bb677c4241463a11d05f29a1f
        ngmix/keep_prepix: 40ba9d5e6afce5422e47c26aaf209cd21e74a5fc
        galsim/#720: 993044f1eb13d577e27d6f48fed23b5059ed8387
    - run-bd37mcal01
        nsim/master: 878af72c52a49b7f12a06f2d68a2a27c43a7d599
        ngmix/keep_prepix: 40ba9d5e6afce5422e47c26aaf209cd21e74a5fc
        galsim/#720: 993044f1eb13d577e27d6f48fed23b5059ed8387
    - run-bd37mcal02
        at bnl
        nsim/master: e94dc404f95fe368f92faf0481e3c9d35cfd06a5
        ngmix/keep_prepix: 8537b9dacde34d75907cb3dc0c92022384470a9f
        galsim/#720: 993044f1eb13d577e27d6f48fed23b5059ed8387
        s/n > 7
            R: [ 0.52606754  0.52637056]
            Rpsf: [ 0.19436534  0.1941443 ]
            Rsel: [ 0.00184014  0.00184935]
            Rpsf_sel: [ 0.00472521  0.00473485]
            without correction
              m1: 3.306e-03 +/- 4.296e-04 c1: -3.796e-06 +/- 1.522e-05  r1: -0.0015
              m2: 3.231e-03 +/- 4.198e-04 c2: 9.404e-05 +/- 1.480e-05  r2: 0.00557
              m:  3.269e-03 +/- 3.006e-04
              c1: -3.794e-06 +/- 1.502e-05
              c2: 9.405e-05 +/- 1.502e-05
            with correction
              m1: -1.923e-04 +/- 4.281e-04 c1: -3.780e-06 +/- 1.516e-05  r1: -0.0015
              m2: -2.814e-04 +/- 4.183e-04 c2: 3.090e-05 +/- 1.475e-05  r2: 0.00557
              m:  -2.367e-04 +/- 2.995e-04
              c1: -3.778e-06 +/- 1.497e-05
              c2: 3.091e-05 +/- 1.497e-05


- sim-stars01
    - stars, goes with bd21

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
    - see aug 19 at slac and aug 20 at bnl
    - the aug 20 run at bnl shows detection of Rpsf in the deep
- run-bd04mcal-t02
    - bnl
    - accidentally ran with prior g sigma 0.3 in both
        the degrade and normal run
    - I had modified it back but did not install it

        meas: 0.0352601 +/- 0.000116722, -0.000132118 +/- 0.000123185
        fracdiff: 7.43e-03 +/- 3.33e-03
    - additive seems to be gone now...

- run-bd04zmcal-degrade03
    - same as run-bd04zmcal-degrade02
    - running at bnl to see if we recover additive
    - short run 100,000
    - I used start_noise_factor of 100, guessing this was
        what I used before.
    - Yes, the psf sens. is recovered
    - why not recovered in some other runs?
        - slac vs bnl?
            - same at slac
        - galsim versions?
            - I'm on master (1.4) at bnl and recovered it with 03
        - start noise?
            - tested with 10 times less start noise, in run-bd04zmcal-degrade04
            - looks similar
        - g prior width 0.3 instead of 0.2
            - run-bd04zmcal-degrade03 also detects and it was 0.2

- run-bd04zmcal-degrade04
    - same as run-bd04zmcal-degrade02 but start noise factor 1000
    - looks similar
- run-bd04zmcal-degrade05
    - same as run-bd04zmcal-degrade03 but run at slac
    - agrees with bnl runs

#  trying integration
- run-bd04zmcal-degint01
    - still additive errors same
    - turns out this was model bias, see below

- sim-bd05
    - high s/n
- run-bd05zmcal-degrade01
    - bnl
- run-bd05mcal-t01
    - bnl
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

- ARGH the galsim versions were different!  I was using master at bnl and 1.3
at slac I'm installing new version of galsim at slac now

- run-ggnr02zmcal-degrade02
    - run at slac with new galsim
- run-ggnr02mcal-t02
    - run at slac with new galsim
    - still additive problems
        meas: 0.035124 +/- 0.000111447, 0.00210291 +/- 0.000111993
        fracdiff: 3.54e-03 +/- 3.18e-03

        mcal_psf_sens 
        [ -5.43e-09, -9.30e-05] +/- [ 7.33e-09, 3.41e-06 ]

        but note the mcal_g from the deep run
            [-0.00019484,  0.00137993] +/- [ 0.00011141,  0.00011147]
        subtracting that gives very small additive (note above it gets
        boosted up by the R division)

- run-ggnr02zmcal-degrade03
    - first of new style sim

- run-ggrn02mcal-t03
    - with deep data from run-ggnr02zmcal-degrade03
        meas: 0.0352673 +/- 0.000113649, -0.0005973 +/- 0.000105899
        fracdiff: 7.64e-03 +/- 3.25e-03
    - raw g mean is 0.0250717, so would need R = 0.716335
        - what to check with short run
            - lowering start noise factor?
            - more nrand?

- run-ggnr02zmcal-degrade04
    - smaller run with nrand=100
    - not enough stats, but looks generally good
    - noticeably less wide than 05

- run-ggnr02zmcal-degrade05
    - smaller run with nrand=40
    - start fac of 1000 instead of 10, still 40 randomizations
    - I was averaging ellipticity with noshear in there too.  Could this
        make a difference at the few in ten thousand level?
    - not enough stats, but looks generally good


- run-ggnr02mcal-t04
    - step 0.02
    - using run-ggnr02zmcal-degrade07
        meas: 0.0350959 +/- 0.000110924, -0.000530868 +/- 0.000110397
        fracdiff: 2.74e-03 +/- 3.17e-03
    - worth longer run
- run-ggnr02zmcal-degrade06
    - another short run just to see if width continues to decrease
    - start-noise-factor 1000
    - nrand 160
    - step 0.02
    - not clear

- run-ggnr02zmcal-degrade07
    - longer version of 06

- run-ggnr02zmcal-degrade-diffn01
    - same_noise: False
    - just add the noise once at the beginning and *then* do metacal,
        rather than adding noise after metacal
    - way too big Rpsf 0.0042
    - way too small R 6.62088972e-01 (need 0.716) 

- run-ggnr02mcal-01
    - up to 40,000,00 from 4,000,000 in -t04
    - error should be 0.001 in fracdev
    - using run-ggnr02zmcal-degrade07 as deep dataset
        meas: 0.0351135 +/- 3.49411e-05, -0.000425756 +/- 3.44229e-05
        fracdiff: 3.24e-03 +/- 9.98e-04

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


pcal
-----
run-gg13pcal-01
    - ran with noise 0.75
    m1: -5.857e-03 +/- 3.623e-03 c1: -1.441e-04 +/- 1.336e-04
    m2: -9.737e-03 +/- 4.353e-03 c2: 5.074e-05 +/- 1.416e-04
run-gg13pcal-02
    - ran with noise 1.0
    m1: -1.876e-02 +/- 3.677e-03 c1: 7.494e-05 +/- 1.356e-04
    m2: -2.961e-02 +/- 4.456e-03 c2: 1.371e-04 +/- 1.449e-04

run-gg13pcal-simp-02
    - ran with 0.85 noise
    - simp correction


ncal
----
see if we can fix noise bias if the model is correct

- run-gg13ncal-01
    - mode s/n ~ 12
    - dumb weights

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




- testing using deep prior sample to simulate the effects of selections
    - sim-bd16 for deep
    - run-bd16zmax-lownoise02 to get priors used to simulate
    - sim-bd16re01 simulated from fits to above
        argh.... noise was too low, s/n distribution is way off
        I should probably just use CM
        80/100
          m1: 1.745e-03 +/- 7.375e-04 c1: -8.016e-06 +/- 2.606e-05  r1: -0.0629
          m2: -4.098e-04 +/- 7.455e-04 c2: -3.932e-05 +/- 2.560e-05  r2: 0.0541
          m:  6.973e-04 +/- 5.252e-04
          c1: -5.687e-06 +/- 2.586e-05
          c2: -3.726e-05 +/- 2.585e-05
        s/n > 20
          557901/1000000
            (436867/1000000 in original.... hmm...)
          m1: 3.055e-04 +/- 7.854e-04 c1: 9.611e-06 +/- 2.775e-05  r1: -0.0629
          m2: -4.309e-04 +/- 8.141e-04 c2: -3.012e-05 +/- 2.796e-05  r2: 0.0541
          m:  -5.252e-05 +/- 5.658e-04
          c1: 1.041e-05 +/- 2.786e-05
          c2: -2.942e-05 +/- 2.785e-05


- sim-bdd26
    - bulge/disk with different ellipticities.  No offset though
    almost all
          m1: 1.409e-03 +/- 6.013e-04 c1: 3.563e-05 +/- 2.047e-05  r1: -0.00251
          m2: 1.909e-03 +/- 6.379e-04 c2: -3.313e-05 +/- 2.219e-05  r2: -0.0197
          m:  1.664e-03 +/- 4.389e-04
          c1: 3.561e-05 +/- 2.136e-05
          c2: -3.296e-05 +/- 2.137e-05

    95/100
          m1: 1.511e-03 +/- 6.208e-04 c1: 3.476e-05 +/- 2.114e-05  r1: -0.00251
          m2: 2.102e-03 +/- 6.524e-04 c2: -3.181e-05 +/- 2.270e-05  r2: -0.0197
          m:  1.813e-03 +/- 4.509e-04
          c1: 3.474e-05 +/- 2.195e-05
          c2: -3.161e-05 +/- 2.195e-05
    T_r > 0.5 psf T_r
            995884/1000000
            sel: [ 1.00018206  0.99996179]
          m1: -7.068e-05 +/- 6.202e-04 c1: 3.362e-05 +/- 2.111e-05  r1: -0.00251
          m2: 1.464e-04 +/- 6.539e-04 c2: -2.522e-05 +/- 2.275e-05  r2: -0.0197
          m:  4.018e-05 +/- 4.512e-04
          c1: 3.361e-05 +/- 2.196e-05
          c2: -2.515e-05 +/- 2.196e-05
    T_r > 0.66 psf T_r
        961200/1000000


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


#
# mike's notes on using a galsim config from within python, rather
# than writing the files
#

# Define the configuration dict
config = {}
config['gal'] = ...
config['psf'] = ...
config['image'] = ... # Just things like noise, wcs, etc. that matter for the postage stamps
# Or can read from a config file 
config = galsim.config.ReadConfig(config_file)

# If necessary, prepare anything that uses an input field:
galsim.config.ProcessInput(config)

# Initialize the rng.  We normally reset the rng for each obj_num, but in this use case, you don't 
# have to.  You can just let it run through the whole simulation if you want.
galsim.config.SetupConfigRNG(config)
# Or if you prefer to set it directly, you can just do this:
config['rng'] = galsim.BaseDeviate(seed)

# Depending on how you want to interact with the image parameters, you might want to 
# run this once here at the top, which will use a common wcs, etc. for all stamps.  Or you might
# want to run it inside the loop with image_num=i, which will re-parse those parameters for
# each stamp.
galsim.config.SetupConfigImageNum(config, image_num=0)

# Similarly, this will define the full image if you want.  But you might want to redo this for each
# stamp with its size.  This is where the wcs actually gets processed, so if you want a different
# wcs for each stamp, this needs to be inside the loop.  Otherwise, fine to put it here.
galsim.config.SetupConfigImageSize(config, full_xsize, full_ysize)

# Now it should all be ready to run through the building of postage stamps:
for i in xrange(ngal):
    # xsize,ysize may be omitted, in which case it will check the config dict for a stamp size 
    # specification, or if none, it will use GalSim's automatic image size for that object.
    stamp, current_var = galsim.config.BuildStamp(config, obj_num=i, xsize=xsize, ysize=ysize)

    # If you want to grab the generated PSF for something, it is stored as
    psf = config['psf']['current_val']
    # Likewise the galaxy profile will be at
    gal = config['gal']['current_val']
