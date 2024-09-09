# Calculates the NLCC contribution to the force
function calc_stress_nlcc!(
    atoms, pspots, pw, xc_calc, xcfunc, Rhoe, rhoe_core, stress_nlcc
)

    Nspecies = atoms.Nspecies
    Npoints = prod(pw.Ns)
    Ng = pw.gvec.Ng
    idx_g2shells = pw.gvec.idx_g2shells
    idx_g2r = pw.gvec.idx_g2r
    G = pw.gvec.G
    G2 = pw.gvec.G2
    G2_shells = pw.gvec.G2_shells
    Ngl = size(G2_shells, 1)

    # XXX: probably add this as field in Hamiltonian or PsPotNL?
    using_nlcc = zeros(Bool, Nspecies)
    for isp in 1:Nspecies
        using_nlcc[isp] = pspots[isp].is_nlcc
    end
    # early return
    if !any(using_nlcc)
        fill!(stress_nlcc, 0.0)
        return
    end  
    #
    # Calculate exchange-correlation potential
    #
    Nspin = 1 # FIXME: need to be deduced from size of rhoe_core
    Vxc = zeros(Float64, Npoints, Nspin)
    epsxc = zeros(Float64, Npoints, Nspin) # not used
    ρ = zeros(Float64, Npoints)
    @views ρ[:] .= Rhoe[:,1] + rhoe_core[:,1]
    if xcfunc == "VWN"
        @views calc_epsxc_Vxc_VWN!( xc_calc, ρ, epsxc[:,1], Vxc[:,1] )
    elseif xcfunc == "PBE"
        @views calc_epsxc_Vxc_PBE!( xc_calc, pw, ρ, epsxc[:,1], Vxc[:,1] )
    else
        error("Unsupported xcfunc")
    end

    ctmp = zeros(ComplexF64, Npoints)
    for ir in 1:Npoints
        ctmp[ir] = Vxc[ir,1]
    end
    #
    planfw = plan_fft!( zeros(ComplexF64,pw.Ns) ) # using default plan
    ff = reshape(ctmp, pw.Ns)
    planfw*ff
    @views ctmp[:] /= Npoints # rescale
    # Now, ctmp is Vxc(G)


    fill!(stress_nlcc, 0.0)
    rhoecgl = zeros(Float64, Ngl)
    drhoecgl = rhoecgl # use rhoecgl memory
    strf = calc_strfact( atoms, pw )
    fact = 1 # FIXME: 2 if using gamma only
    #
    sigmadiag = 0.0
    #
    for isp in 1:Nspecies
        #
        psp = pspots[isp]        
        # skip if no nlcc
        if !psp.is_nlcc
            continue
        end
        #
        PWDFT._calc_rhoecgl!(psp, G2_shells, rhoecgl)
        # Remember that rhoecgl does not include 1/CellVolume factor
        #
        # G = 0 contribution
        sigmadiag += real( conj(ctmp[1])*strf[1,isp] ) * rhoecgl[1] / pw.CellVolume
        # G != 0 contribution
        for ig in 2:Ng
            igl = idx_g2shells[ig]
            ip = idx_g2r[ig]
            crg = real( conj(ctmp[ip]) * strf[ig,isp] )
            sigmadiag += crg * rhoecgl[igl] * fact / pw.CellVolume
        end
        #
        _calc_drhoecgl!(psp, G2_shells, drhoecgl)
        # Non diagonal term
        for ig in 2:Ng, l in 1:3, m in 1:3
            igl = idx_g2shells[ig]
            ip = idx_g2r[ig]
            crg = real( conj(ctmp[ip]) * strf[ig,isp] )
            GlGm = G[l,ig] * G[m,ig] * fact / sqrt(G2[ig])
            stress_nlcc[l,m] += crg * drhoecgl[igl] * GlGm / pw.CellVolume
        end
    end
    # Add diagonal component
    for l in 1:3
        stress_nlcc[l,l] += sigmadiag
    end
    return
end



function _calc_drhoecgl!(
    psp::PsPot_UPF,
    G2_shells::Vector{Float64},
    drhoecgl::Vector{Float64}
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

    Ngl = length(G2_shells)
    aux = zeros(Float64, Nr)
    pref = 4π
    # G = 0 term
    if G2_shells[1] < 1e-8
        drhoecgl[1] = 0.0
    end
    #
    # G != 0 term
    for igl in 2:Ngl
        Gx = sqrt(G2_shells[igl])
        for ir in 1:Nr
            aux[ir] = r[ir] * rho_atc[ir]*( r[ir] * cos(Gx*r[ir])/Gx - sin(Gx*r[ir])/Gx^2 )
        end
        drhoecgl[igl] = pref*PWDFT.integ_simpson(Nr, aux, rab)
    end
    return
end