function calc_rhoe_core!(
    atoms::Atoms,
    pw::PWGrid,
    pspots::Vector{PsPot_UPF},
    rhoe_core
)
    # XXX: rhoe_core is assumed to have shape (Npoints,1)

    Nspecies = atoms.Nspecies
    Npoints = prod(pw.Ns)
    G2_shells = pw.gvec.G2_shells
    Ngl = length(G2_shells)
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r
    idx_g2shells = pw.gvec.idx_g2shells
    CellVolume = pw.CellVolume

    Nspin = 1 # FIXME
    strf = calc_strfact(atoms, pw)
    rhoecgl = zeros(Float64, Ngl)
    rhoecG = zeros(ComplexF64, Npoints)

    for isp in 1:Nspecies
        psp = pspots[isp]
        # Skip this pspot if nlcc is not included
        if !psp.is_nlcc
            continue
        end
        _calc_rhoecgl!(psp, G2_shells, rhoecgl)
        #println("sum of rhoecgl = ", sum(rhoecgl)/CellVolume)
        for ig in 1:Ng
            ip = idx_g2r[ig]
            igl = idx_g2shells[ig]
            # accumulate
            rhoecG[ip] += strf[ig,isp] * rhoecgl[igl] / CellVolume
        end
    end

    G_to_R!(pw, rhoecG)
    for ip in 1:Npoints
        rhoe_core[ip,1] = real(rhoecG[ip]) * Npoints
    end

    rhoneg = 0.0
    for ip in 1:Npoints
        rhoneg += min(0.0, rhoe_core[ip])
        rhoe_core[ip] = max(0.0, rhoe_core[ip])
    end
    println("Check negative rhoe_core (mean) = ", rhoneg/Npoints)
    println("integ rhoe_core = ", sum(rhoe_core)*CellVolume/Npoints)
    return
end


function _calc_rhoecgl!(
    psp::PsPot_UPF,
    G2_shells::Vector{Float64},
    rhoecgl::Vector{Float64}
)
    r = psp.r
    rab = psp.rab
    rho_atc = psp.rho_atc

    Nr_full = psp.Nr
    Nr = Nr_full
    RCUT = 10.0
    for i in 1:Nr_full
        if r[i] > RCUT
            Nr = i
            break
        end 
    end
    # Force Nr to be odd number
    Nr = 2*floor(Int64, (Nr + 1)/2) - 1
    println("Nr = ", Nr)

    Ngl = length(G2_shells)
    aux = zeros(Float64, Nr)

    pref = 4π

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
            #aux[ir] = r[ir]^2 * rho_atc[ir] * sphericalbesselj(0, Gx*r[ir])
            # or
            aux[ir] = r[ir] * rho_atc[ir] * sin(Gx*r[ir])/Gx
        end
        rhoecgl[igl] = pref*integ_simpson(Nr, aux, rab)
    end
    return
end