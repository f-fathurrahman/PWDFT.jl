# Calculates the NLCC contribution to the force
function calc_forces_nlcc!( Ham::Hamiltonian{PsPot_UPF}, F_nlcc )
    calc_forces_nlcc!(
        Ham.atoms, Ham.pspots, Ham.pw,
        Ham.xc_calc, Ham.xcfunc, Ham.rhoe, Ham.rhoe_core,
        F_nlcc
    )
    return
end


function calc_forces_nlcc!(
    atoms::Atoms,
    pspots::Vector{PsPot_UPF},
    pw, xc_calc, xcfunc, Rhoe, rhoe_core, F_nlcc
)

    atm2species = atoms.atm2species
    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    X = atoms.positions
    Npoints = prod(pw.Ns)
    Ng = pw.gvec.Ng
    idx_g2shells = pw.gvec.idx_g2shells
    idx_g2r = pw.gvec.idx_g2r
    G = pw.gvec.G
    G2_shells = pw.gvec.G2_shells
    Ngl = size(G2_shells, 1)

    # XXX: probably add this as field in Hamiltonian or PsPotNL?
    using_nlcc = zeros(Bool, Nspecies)
    for isp in 1:Nspecies
        using_nlcc[isp] = pspots[isp].is_nlcc
    end
    # early return
    if !any(using_nlcc)
        fill!(F_nlcc, 0.0)
        return
    end

    fact = 1 # FIXME: 2 if using gamma only
  
    #
    # Calculate exchange-correlation potential
    #
    Nspin = 1 # FIXME: need to be deduced from size of rhoe_core
    Vxc = zeros(Float64, Npoints, Nspin)
    epsxc = zeros(Float64, Npoints, Nspin) # not used
    println("sum rhoe_core = ", sum(rhoe_core))
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
    println("sum ctmp before forward FFT: ", sum(ctmp))

    # This will will segfault, FFT plan is not initialized because we read Ham
    # from Julia serialization
    R_to_G!(pw, ctmp)
    #
    # We use resort to the usual way
    #planfw = plan_fft!( zeros(ComplexF64,pw.Ns) ) # using default plan
    #ff = reshape(ctmp, pw.Ns)
    #planfw*ff
    @views ctmp[:] /= Npoints # rescale
    # Now, ctmp is Vxc(G)

    println("sum ctmp after forward FFT: ", sum(ctmp))

    fill!(F_nlcc, 0.0)
    rhoecgl = zeros(Float64, Ngl)
    # core correction term: sum on g of omega*ig*exp(-i*r_i*g)*n_core(g)*vxc
    # G = 0 term gives no contribution
    
    for isp in 1:Nspecies

        psp = pspots[isp]        
        # skip if no nlcc
        if !psp.is_nlcc
            continue
        end
        #
        PWDFT._calc_rhoecgl!(psp, G2_shells, rhoecgl)
        println("sum rhoecgl/pw.CellVolume = ", sum(rhoecgl)/pw.CellVolume)
        # rhoecgl does not include 1/CellVolume factor
        #        
        for ia in 1:Natoms
            if atm2species[ia] != isp
                continue
            end
            for ig in 2:Ng
                igl = idx_g2shells[ig]
                ip = idx_g2r[ig]
                GX = G[1,ig]*X[1,ia] + G[2,ig]*X[2,ia] + G[3,ig]*X[3,ia]
                Sf = sin(GX) + im*cos(GX)
                F_nlcc[:,ia] .+= fact * G[:,ig] * rhoecgl[igl] * real(conj(ctmp[ip]) * Sf)
                # pw.CellVolume factor should cancel out
            end
        end
    end
    return
end
