import numpy

def deblend(mb_obslist):
    """
    assume central is first, deblend the others
    """
    import deblender
    
    nobj = len(mb_obslist)

    # only need the one image, the others are duplicates
    obslist_cen = mb_obslist[0]
    obs_cen = obslist_cen[0]

    image = obs_cen.image
    weight = obs_cen.weight

    cens = [ obsl[0].jacobian.get_cen() for obsl in mb_obslist]

    cens = [ (c[1],c[0]) for c in cens]

    rimage = image[numpy.newaxis,:,:]
    rweight = weight[numpy.newaxis, :, :]

    result = deblender.nmf.deblend(
        rimage,
        weights=rweight,
        peaks=cens,
        strict_constraints="M",
        constraints="MS",
        max_iter=500,
        e_rel=[1e-6,1e-3],
        psf_thresh=3e-3,
        l0_thresh=.005,
        algorithm="GLMM",
        traceback=False,
        convergence_func=None,
        als_max_iter=20,
        monotonicUseNearest=False,
    )
    A, S, model, P_, Tx, Ty, errors = result

    deblended_image = image.copy()
    
    for k in xrange(nobj):
        res = deblender.nmf.get_peak_model(
            A[:,k],
            S[k].flatten(),
            Tx[k],
            Ty[k],
            shape=(S[k].shape),
        )

        tmodel = res[0]

        deblended_image -= tmodel

    obs_cen.image_orig = image
    obs_cen.set_image(deblended_image)

    return obslist_cen

