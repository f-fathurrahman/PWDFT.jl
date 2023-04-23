# Use the density produced by sum_rad_rho to compute xc potential
# and energy, as xc functional is not diagonal on angular momentum
# numerical integration is performed.
function PAW_xc_potential!(
    AE::Bool, ia,
    atoms, pspots, pspotNL,
    xc_calc,
    rho_lm, v_lm
)
  
    isp = atoms.atm2species[ia]
    r2 = pspots[isp].r.^2
    if AE
        rho_core = pspots[isp].paw_data.ae_rho_atc
    else
        rho_core = pspots[isp].rho_atc
    end
    nx = pspotNL.paw.spheres[isp].nx

    Nrmesh = pspots[isp].Nr
    Nspin = size(rho_lm, 3)

    # This will hold the "true" charge density, without r**2 or other factors
    rho_loc = zeros(Float64, Nrmesh, Nspin)

    # radial potential (to be integrated)
    v_rad = zeros(Float64, Nrmesh, nx, Nspin)

    #
    rho_rad = zeros(Float64, Nrmesh, Nspin) 
    arho = zeros(Float64, Nrmesh, 2) # XXX: Nspin = 2
    # it is also 2 in case of noncollinear spin (?)
    
    energy = 0.0
    e_rad = zeros(Float64, Nrmesh)

    for ix in 1:nx

        # LDA (and LSDA) part (no gradient correction)
        # convert _lm density to real density along ix
        PAW_lm2rad!( ia, ix, atoms, pspots, pspotNL, rho_lm, rho_rad )

        # compute the potential along ix
        if Nspin == 2
            for k in 1:Nrmesh
                rho_loc[k,1] = rho_rad[k,1]/r2[k]
                rho_loc[k,2] = rho_rad[k,2]/r2[k]
            end
        else
            for k in 1:Nrmesh
                rho_loc[k,1] = rho_rad[k,1]/r2[k]
            end
        end

        #
        # Integrate to obtain the energy
        #
        if Nspin == 1
            @views arho[:,1] .= rho_loc[:,1] .+ rho_core[:]
            @views calc_epsxc_Vxc_VWN!( xc_calc, arho[:,1], e_rad[:,1], v_rad[:,ix,1] )
            e_rad .= e_rad .* ( rho_rad[:,1] .+ rho_core .* r2 )
        else
            # This is not yet working
            println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            println("WARNING: This is not yet tested")
            println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            @views arho[:,1] .= rho_loc[:,1] .+ rho_loc[:,2] .+ rho_core[:]
            @views arho[:,2] .= rho_loc[:,1] .- rho_loc[:,2]
            calc_epsxc_Vxc_VWN!( xc_calc, arho, e_rad[:,:], v_rad[:,ix,:] )
            e_rad .= e_rad .* ( rho_rad[:,1] .+ rho_rad[:,2] .+ rho_core .* r2 )
        end
    
        # Integrate to obtain the energy
        wx = pspotNL.paw.spheres[isp].ww[ix]
        energy += wx*PWDFT.integ_simpson( Nrmesh, e_rad, pspots[isp].rab )
  
    end

    #
    # Recompose the sph. harm. expansion
    lmax_loc = pspots[isp].lmax_rho + 1
    PAW_rad2lm!( ia, atoms, pspotNL, lmax_loc, v_rad, v_lm )

    #
    # Add gradient correction, if necessary
    #
    # Not yet implemented

    return energy

end