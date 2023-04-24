function update_from_rhoe!(Ham, psiks::BlochWavefunc, Rhoe)
    # RhoeG is not given we need to calculate it
    RhoeG = _rhoeG_from_rhoe(Ham, Rhoe)
    #
    return update_from_rhoe!(Ham, psiks, Rhoe, RhoeG)
end

function _rhoeG_from_rhoe(Ham, Rhoe)
    Nspin = size(Rhoe,2)
    Npoints = size(Rhoe,1)
    RhoeG = zeros(ComplexF64, Npoints, Nspin)
    ctmp = zeros(ComplexF64, Npoints)
    for ispin in 1:Nspin
        ctmp[:] = Rhoe[:,1]
        #
        R_to_G!(Ham.pw, ctmp) # FIXME: add method
        RhoeG[:,1] = ctmp[:]/Npoints
        # Need to be careful about the normalization
        # Check the charge in G-space
        charge = RhoeG[1,1]*Ham.pw.CellVolume
        # println("Check charge from RhoeG: ", charge)    
    end
    return RhoeG
end


# Can be used for metaGGA (if psiks is not nothing)
# FIXME: Think of better API for this
# Probably by setting psiks as optional argument
function update_from_rhoe!(Ham, psiks, Rhoe, RhoeG)

    Ham.rhoe[:,:] = Rhoe[:,:] # Need copy?

    # Save old potential
    Ham.potentials.TotalOld[:,:] .= Ham.potentials.Total[:,:]

    # Reset total effective potential to zero
    fill!(Ham.potentials.Total, 0.0)

    # FIXME: also set Ham.potentials.XC and Ham.potentials.Hartree
    Exc, Evtxc = _add_V_xc!( Ham, psiks, Rhoe, RhoeG )
    Ehartree = _add_V_Hartree!( Ham, Rhoe, RhoeG )

    # Add V_Ps_loc contribution
    Nspin = Ham.electrons.Nspin
    for ispin in 1:Nspin
        Ham.potentials.Total[:,ispin] .+= Ham.potentials.Ps_loc[:]
    end

    if Ham.pw.using_dual_grid
        dense_to_smooth!( Ham.pw, Ham.potentials.Total, Ham.potentials.TotalSmooth )
    end
    # XXX: Also need to interpolate xc_calc.Vtau or kedtau in case of USPP

    # PAW stuffs
    if any(Ham.pspotNL.are_paw)
        # becsum is assumed to be calculated elsewhere
        # Possibly alongside Rhoe calculation
        EHxc_paw = PAW_potential!( Ham )
        # EHxc_paw is currently not returned
        # We set is here
        Ham.pspotNL.paw.EHxc_paw = EHxc_paw
    end

    # Also update nonlocal potential coefficients here
    calc_newDeeq!( Ham )
    # PAW-specific stuffs is also updated here


    return Ehartree, Exc, Evtxc # energies?
end


# Can be used for metaGGA functionals
# FIXME: allow psiks to be nothing (?)
function _add_V_xc!(Ham, psiks, Rhoe, RhoeG)

    Nspin = Ham.electrons.Nspin
    
    @assert Nspin == 1

    Npoints = size(Rhoe, 1)
    epsxc = zeros(Float64, Npoints)
    Vxc = Ham.potentials.XC

    dVol = Ham.pw.CellVolume / Npoints

    # XC potential
    if Ham.rhoe_core == nothing
        #
        # No core-correction
        #
        if Ham.xcfunc == "SCAN"

            @views calc_epsxc_Vxc_SCAN!( Ham, psiks, Rhoe[:,1], epsxc, Vxc[:,1] )
            # Note that we are using the inplace version of calc_epsxc_Vxc_SCAN
        
        elseif Ham.xcfunc == "PBE"
        
            @views calc_epsxc_Vxc_PBE!( Ham, Rhoe[:,1], epsxc, Vxc[:,1] )
        
        else
            # VWN
            epsxc[:], Vxc[:,1] = calc_epsxc_Vxc_VWN( Ham.xc_calc, Rhoe[:,1] )
        end
        Exc = sum(epsxc .* Rhoe[:,1])*dVol
    else
        #
        # Using core-correction
        #
        if Ham.xcfunc == "VWN"
            epsxc[:], Vxc[:,1] = calc_epsxc_Vxc_VWN( Ham.xc_calc, Rhoe[:,1] + Ham.rhoe_core )
            Exc = sum(epsxc .* (Rhoe[:,1] + Ham.rhoe_core))*dVol
        else
            println("Other than VWN, core correction is yet not supported")
            error()
        end
    end

    # Also calculate Evtxc
    Evtxc = sum(Vxc[:,1] .* Rhoe)*dVol # Evtxc does not include rhoe_core

    # XXX: Evtxc is vtxc in QE, it seems that is is not used for total energy calculation

    Ham.potentials.Total[:,1] += Vxc[:,1] # Update

    return Exc, Evtxc
end


# Note that RhoeG is already in FFT grid
function _add_V_Hartree!(Ham, Rhoe, RhoeG)

    pw = Ham.pw
    gvec = pw.gvec

    G2 = gvec.G2
    Ng = gvec.Ng
    idx_g2r = gvec.idx_g2r
    Npoints = prod(pw.Ns) # dense

    VhG = zeros(ComplexF64, Npoints)
    Ehartree = 0.0
    # skip ig = 1, it is set to zero
    for ig in 2:Ng
        ip = idx_g2r[ig]
        Ehartree = Ehartree + 2π*( real(RhoeG[ip])^2 + imag(RhoeG[ip])^2 )/G2[ig]
        VhG[ip] = 4π * RhoeG[ip]/G2[ig]
    end

    Ehartree *= Ham.pw.CellVolume

    G_to_R!(pw, VhG)
    VhG[:] *= Npoints # XXX: scale by Npoints
    Ham.potentials.Hartree[:] = real(VhG) # update
    Ham.potentials.Total[:,1] += real(VhG)

    return Ehartree
end