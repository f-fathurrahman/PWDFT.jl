function calc_forces_nlcc( Ham::Hamiltonian{PsPot_UPF} )
    F_nlcc = zeros(Float64, 3, Ham.atoms.Natoms)
    calc_forces_nlcc!(Ham, F_nlcc)
    return F_nlcc
end

# Calculates the NLCC contribution to the force
function calc_forces_nlcc!( Ham::Hamiltonian{PsPot_UPF}, F_nlcc )
    calc_forces_nlcc!(
        Ham.atoms, Ham.pspots, Ham.pw,
        Ham.xc_calc, Ham.xcfunc, Ham.rhoe, Ham.rhoe_core,
        F_nlcc
    )
    return
end


# For PsPot_GTH, the contributions are zeros. We do not yet
# support nonlinear core corrections in PsPot_GTH
function calc_forces_nlcc!( Ham::Hamiltonian{PsPot_GTH}, F_nlcc )
    fill!(F_nlcc, 0.0)
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
    Nspin = size(rhoe_core, 2)
    @assert Nspin == 1
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

    # Bring to G-space
    R_to_G!(pw, ctmp)
    @views ctmp[:] /= Npoints # rescale
    # Now, ctmp is Vxc(G)

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
        _calc_rhoecgl!(psp, G2_shells, rhoecgl)
        # Reminder: rhoecgl does not include 1/CellVolume factor
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
