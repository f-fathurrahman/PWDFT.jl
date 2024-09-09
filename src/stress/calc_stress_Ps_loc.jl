function calc_stress_Ps_loc!( atoms, pw, pspots, Rhoe, stress_Ps_loc )

    fill!(stress_Ps_loc, 0.0)

    Nspecies = atoms.Nspecies

    Npoints = prod(pw.Ns)
    G = pw.gvec.G
    Ng = pw.gvec.Ng
    G2_shells = pw.gvec.G2_shells
    Ngl = length(G2_shells)
    idx_g2r = pw.gvec.idx_g2r
    idx_g2shells = pw.gvec.idx_g2shells
    CellVolume = pw.CellVolume

    ctmp = zeros(ComplexF64, Npoints)
    # Use Nspin=1
    for ip in 1:Npoints
        ctmp[ip] = Rhoe[ip,1]
    end
    # We create again plan_fft
    planfw = plan_fft!(zeros(ComplexF64, pw.Ns)) # using default plan
    ff = reshape(ctmp, pw.Ns)
    planfw*ff
    @views ctmp[:] /= Npoints # rescale

    fact = 1.0 # it is 2 if using gamma only

    Vgl = zeros(Float64, Ngl)
    strf = calc_strfact( atoms, pw )
    evloc = 0.0
    for isp in 1:Nspecies
        eval_Vloc_G!( pspots[isp], G2_shells, Vgl )
        # G = 0 contribution
        evloc += real(ctmp[1] * strf[1,isp]) * Vgl[1] / CellVolume
        # G != 0 contribution
        for ig in 2:Ng
            ip = idx_g2r[ig]
            igl = idx_g2shells[ig]
            crg = real(conj(ctmp[ip])*strf[ig,isp])
            evloc += crg * Vgl[igl] * fact / CellVolume
        end
    end
    # skip 2d
    println("evloc = ", evloc)

    dVgl = zeros(Float64, Ngl)
    for isp in 1:Nspecies
        eval_dVloc_G!( pspots[isp], G2_shells, dVgl )
        # no G=0 contribution
        for ig in 2:Ng, l in 1:3, m in 1:l
            ip = idx_g2r[ig]
            igl = idx_g2shells[ig]
            GlGm = G[l,ig] * G[m,ig] * fact
            crg = real(conj(ctmp[ip])*strf[ig,isp])
            stress_Ps_loc[l,m] += crg * 2.0 * dVgl[igl] * GlGm / CellVolume
        end
    end
    #
    for l in 1:3
        stress_Ps_loc[l,l] += evloc
        for m = 1:(l-1)
            stress_Ps_loc[m,l] = stress_Ps_loc[l,m]
        end
    end
    #
    return
end
