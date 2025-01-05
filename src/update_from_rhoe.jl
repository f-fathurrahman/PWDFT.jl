function update_from_rhoe!(Ham, psiks::BlochWavefunc, Rhoe)
    # RhoeG is not given, we need to calculate it first
    RhoeG = _rhoeG_from_rhoe(Ham, Rhoe)
    #
    return update_from_rhoe!(Ham, psiks, Rhoe, RhoeG)
end

# Calculate reciprocal space representation of Rhoe
# This simply calls R_to_G!
function _rhoeG_from_rhoe(Ham, Rhoe)
    Nspin = size(Rhoe,2)
    Npoints = size(Rhoe,1)
    RhoeG = zeros(ComplexF64, Npoints, Nspin)
    ctmp = zeros(ComplexF64, Npoints)
    for ispin in 1:Nspin
        ctmp[:] .= Rhoe[:,ispin]
        #
        R_to_G!(Ham.pw, ctmp) # FIXME: add method
        RhoeG[:,ispin] = ctmp[:]/Npoints
    end
    # Need to be careful about the normalization
    # Check the charge in G-space
    #charge = RhoeG[1,1]*Ham.pw.CellVolume
    # println("Check charge from RhoeG: ", charge)    
    return RhoeG
end


# Can be used for metaGGA (if psiks is not nothing)
# FIXME: Think of better API for this
# Probably by setting psiks as optional argument
function update_from_rhoe!(Ham, psiks, Rhoe, RhoeG)

    Nspin = Ham.electrons.Nspin
    Npoints = size(Rhoe, 1)

    # Save old potential
    # This will be used to calculate SCF correction to the forces
    #Ham.potentials.TotalOld[:,:] .= Ham.potentials.Total[:,:]
    # XXX: We might only need to save this at the end of a
    # converged SCF calculation, or, calculate it as needed

    # Save old potential (only Hartree and XC)
    for ispin in 1:Nspin, ip in 1:Npoints
        Ham.potentials.TotalOld[ip,ispin] = Ham.potentials.XC[ip,ispin] + Ham.potentials.Hartree[ip]
    end
    # This might be useful to potential mixing, too.
    # XXX: Rename TotalOld to HxcPrev ?

    # Reset total effective potential to zero
    fill!(Ham.potentials.Total, 0.0)

    # Update Ham.potentials.XC and Ham.potentials.Hartree
    # and also compute Exc and Ehartree
    Exc = _add_V_xc!( Ham, psiks, Rhoe )
    Ehartree = _add_V_Hartree!( Ham, RhoeG )

    # Add V_Ps_loc contribution
    for ispin in 1:Nspin
        Ham.potentials.Total[:,ispin] .+= Ham.potentials.Ps_loc[:]
    end

    if Ham.pw.using_dual_grid
        for ispin in 1:Nspin
            @views dense_to_smooth!( Ham.pw, Ham.potentials.Total[:,ispin], Ham.potentials.TotalSmooth[:,ispin] )
        end
    end
    # XXX: Also need to interpolate xc_calc.Vtau or kedtau in case of USPP

    # PAW stuffs
    if any(Ham.pspotNL.are_paw)
        # becsum is assumed to be calculated elsewhere
        # Possibly alongside Rhoe calculation
        EHxc_paw = PAW_potential!( Ham )
        #println("EHxc_paw = ", EHxc_paw)
        #println("sum ddd_paw = ", sum(Ham.pspotNL.paw.ddd_paw))
        #println("sum becsum = ", sum(Ham.pspotNL.becsum))
        # EHxc_paw is currently not returned
        # We set is here
        Ham.pspotNL.paw.EHxc_paw = EHxc_paw
        #
        PAW_symmetrize_ddd!( Ham )
    end

    # Also update nonlocal potential coefficients here
    calc_newDeeq!( Ham )
    # PAW-specific stuffs is also updated here


    return Ehartree, Exc
end


# Can be used for metaGGA functionals
# FIXME: allow psiks to be nothing (?)
function _add_V_xc!(Ham, psiks, Rhoe)

    Nspin = Ham.electrons.Nspin

    Npoints = size(Rhoe, 1)
    epsxc = zeros(Float64, Npoints)
    Vxc = Ham.potentials.XC # This will be written

    dVol = Ham.pw.CellVolume / Npoints

    # XC potential
    # XXX: This can be simplified now that Ham.rhoe_core can be zeros
    if isnothing(Ham.rhoe_core)
        #
        # No core-correction
        #
        if Ham.xcfunc == "SCAN"
            
            @assert Nspin == 1
            @views calc_epsxc_Vxc_SCAN!( Ham, psiks, Rhoe[:,1], epsxc, Vxc[:,1] )
            # Note that we are using the inplace version of calc_epsxc_Vxc_SCAN
        
        elseif Ham.xcfunc == "PBE"

            @assert Nspin == 1
            @views calc_epsxc_Vxc_PBE!( Ham.xc_calc, Ham.pw, Rhoe[:,1], epsxc, Vxc[:,1] )
        
        else
            # VWN
            calc_epsxc_Vxc_VWN!( Ham.xc_calc, Rhoe, epsxc, Vxc )
        end
        Exc = sum(epsxc .* Rhoe)*dVol # is this working for spinpol?
    
    else

        #
        # Using core-correction
        #

        if Nspin == 2
            Rhoe[:,1] .+= Ham.rhoe_core*0.5
            Rhoe[:,2] .+= Ham.rhoe_core*0.5
        else
            Rhoe[:,1] .+= Ham.rhoe_core
        end

        if Ham.xcfunc == "VWN"
            calc_epsxc_Vxc_VWN!( Ham.xc_calc, Rhoe, epsxc, Vxc )
        elseif Ham.xcfunc == "PBE"
            # check if this is working
            calc_epsxc_Vxc_PBE!( Ham.xc_calc, Ham.pw, Rhoe, epsxc, Vxc )
        else
            # This is SCAN
            error("Core correction is yet not supported in SCAN")
        end
        Exc = sum(epsxc .* Rhoe)*dVol

        # Recover
        if Nspin == 2
            Rhoe[:,1] .-= Ham.rhoe_core*0.5
            Rhoe[:,2] .-= Ham.rhoe_core*0.5
        else
            Rhoe[:,1] .-= Ham.rhoe_core
        end

    end

    # Also calculate Evtxc
    #Evtxc = sum(Vxc[:,1] .* Rhoe)*dVol # Evtxc does not include rhoe_core
    # term in GGA is not accounted here (ref gradcorr.f90)
    #
    # vtxcgc = vtxcgc - SUM( dh(:) * rhoaux(:,is) )
    #
    # Anyway, Evtxc is not used in our calculation

    # XXX: Evtxc is vtxc in QE, it seems that it is not used for total energy calculation
    # Evtxc will be used for stress calculation

    Ham.potentials.Total .+= Vxc # Update

    return Exc
end


# Note that RhoeG is already in FFT grid
function _add_V_Hartree!(Ham, RhoeG)

    pw = Ham.pw
    gvec = pw.gvec

    G2 = gvec.G2
    Ng = gvec.Ng
    idx_g2r = gvec.idx_g2r
    Npoints = prod(pw.Ns) # dense
    Nspin = size(RhoeG, 2)
    # get total RhoeG
    if Nspin == 2
        ρG = RhoeG[:,1] + RhoeG[:,2]
    else
        ρG = RhoeG[:,1]
    end
    VhG = zeros(ComplexF64, Npoints)
    
    Ehartree = 0.0
    # skip ig = 1, it is set to zero
    for ig in 2:Ng
        ip = idx_g2r[ig]
        # compute Ehartree
        Ehartree += 2π*( real(ρG[ip])^2 + imag(ρG[ip])^2 )/G2[ig]
        # and Hartree potential in G-space
        VhG[ip] = 4π * ρG[ip]/G2[ig]
    end

    # Transform Hartree potential to R-space
    G_to_R!(pw, VhG)
    VhG[:] *= Npoints # XXX: scale by Npoints
    Ham.potentials.Hartree[:] = real(VhG) # update
    Ham.potentials.Total .+= real(VhG)

    # Hartree energy will be returned
    Ehartree *= Ham.pw.CellVolume
    return Ehartree
end
