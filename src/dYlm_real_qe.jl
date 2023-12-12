function dYlm_real_qe!(
    lmax::Int64,
    R::Matrix{Float64},
    dYlm::Matrix{Float64},
    ipol::Int64
)

    # dYlm is the output
    # ipol is direction, input (x=1,y=1,z=1)

    ngy = size(R, 2)
    nylm = size(dYlm, 2)

    fill!(dYlm, 0.0) # zero out the input

    # nylm, ngy, g, gg, dylm, ipol)
    # compute ∂ Y_lm(G) / ∂ (G)_ipol
    # using simple numerical derivation (SdG)
    # The spherical harmonics are calculated in ylmr2

    Δ = 1e-6
    # dg is the finite increment for numerical derivation:
    # dg = delta |G| = delta * sqrt(gg)
    # dgi= 1 /(delta * sqrt(gg))
    # gx = g +/- dg
    # ggx = gx^2
    #
    gx = zeros(Float64, 3, ngy)
    ggx = zeros(Float64, 3, ngy)
    dg = zeros(Float64, ngy)
    dgi = zeros(Float64, ngy)
    ylmaux = zeros(Float64, ngy, nylm)

    for ig in 1:ngy
        gg = R[1,ig]^2 + R[2,ig]^2 + R[3,ig]^2
        dg[ig] = Δ * sqrt(gg)
        if gg > 1.e-9
            dgi[ig] = 1.0 / dg[ig]
        else
            dgi[ig] = 0.0
        end
    end

    gx[:,:] .= R[:,:] # XXX Need this???
    #call dcopy (3 * ngy, g, 1, gx, 1)

    for ig in 1:ngy
        gx[ipol, ig] = R[ipol,ig] + dg[ig]
    end
    Ylm_real_qe!(lmax, gx, dYlm)

    for ig in 1:ngy
        gx[ipol,ig] = R[ipol,ig] - dg[ig]
    end
    Ylm_real_qe!(lmax, gx, ylmaux)

    # Compute the difference
    @views dYlm[:,:] .= dYlm[:,:] .- ylmaux[:,:]
    #call daxpy(ngy * nylm, - 1.d0, ylmaux, 1, dylm, 1)

    # Using centered difference, multiply dy factor (1/(2*dg))
    for lm in 1:nylm
        for ig in 1:ngy
            dYlm[ig,lm] *= 0.5 * dgi[ig]
        end
    end

    return
end
