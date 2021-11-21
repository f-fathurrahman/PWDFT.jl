function calc_rhoe_core!(atoms, pw, pspots, rhoe_core)

    Nspecies = atoms.Nspecies
    Npoints = prod(pw.Ns)
    G2_shells = pw.gvec.G2_shells
    Ngl = length(G2_shells)
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r
    idx_g2shells = pw.gvec.idx_g2shells
    CellVolume = pw.CellVolume

    rhoe_core = zeros(Float64, Npoints)
    strf = calc_strfact(atoms, pw)
    rhoecgl = zeros(Float64, Ngl)
    rhoecG = zeros(ComplexF64, Npoints)

    neg_rhoec = 0.0
    for isp in 1:Nspecies
        psp = pspots[isp]
        if !psp.is_nlcc
            continue
        end
        _calc_rhoecgl!(psp, G2_shells, rhoecgl)
        for ig = 1:Ng
            ip = idx_g2r[ig]
            igl = idx_g2shells[ig]
            rhoecG[ip] = strf[ig,isp] * rhoecgl[igl]
        end
        #
        rhoe_core[:] = rhoe_core[:] + real(G_to_R(pw, rhoecG)) * Npoints / CellVolume
        for ip in 1:Npoints
            if rhoe_core[ip] < 0.0
                neg_rhoec = neg_rhoec + rhoe_core[ip]
            end
        end
    end
    println("integ rhoec = ", sum(rhoe_core)*CellVolume/Npoints)
    println("neg_rhoec = ", neg_rhoec*CellVolume/Npoints)
    return rhoe_core
end


function _calc_rhoecgl!(
    psp::PsPot_UPF,
    G2_shells::Vector{Float64},
    rhoecgl::Vector{Float64}
)
    r = psp.r
    rab = psp.rab
    Nr = psp.Nr
    rho_atc = psp.rho_atc

    Ngl = length(G2_shells)
    aux = zeros(Float64, Nr)

    pref = 4*pi

    # G=0 term
    if G2_shells[1] < 1e-8
        for ir in 1:Nr
            aux[ir] = r[ir]^2 * rho_atc[ir]
        end
        rhoecgl[1] = pref*integ_simpson(Nr, aux, rab)
    end

    # G != 0 term
    for igl in 2:Ngl
        Gx = sqrt(G2_shells[igl])
        for ir in 1:Nr
            aux[ir] = r[ir]^2 * rho_atc[ir] * sin(Gx*r[ir])/Gx
        end
        rhoecgl[igl] = pref*integ_simpson(Nr, aux, rab)
    end
    return
end