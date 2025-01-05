# Use the density produced by sum_rad_rho to compute xc potential
# and energy, as xc functional is not diagonal on angular momentum
# numerical integration is performed.
function PAW_xc_potential_LDA!(
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
    #v_rad = zeros(Float64, Nrmesh, nx, Nspin)
    v_rad = Vector{Matrix{Float64}}(undef, nx)
    for ix in 1:nx
        v_rad[ix] = zeros(Float64, Nrmesh, Nspin)
    end

    #
    rho_rad = zeros(Float64, Nrmesh, Nspin) 
    arho = zeros(Float64, Nrmesh, Nspin) # XXX: Nspin = 2
    # it is also 2 in case of noncollinear spin (?)
    
    energy = 0.0
    e_rad = zeros(Float64, Nrmesh)

    for ix in 1:nx

        #println()
        #println("PAW_xc_potential: ix = ", ix)

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

        #println("sum rho_loc = ", sum(rho_loc))
        #println("sum rho_core =  ", sum(rho_core))

        #
        # Integrate to obtain the energy
        #
        if Nspin == 1
            @views arho[:,1] .= rho_loc[:,1] .+ rho_core[:]
            # rho_loc excludes the 1/r^2 factor
            #println("sum arho = ", sum(arho))
            @views calc_epsxc_Vxc_VWN!( xc_calc, arho[:,1], e_rad[:,1], v_rad[ix][:,1] )
            @views e_rad .= e_rad .* ( rho_rad[:,1] .+ rho_core .* r2 )
        else
            # This is not yet working
            #println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            #println("WARNING: This is not yet tested")
            #println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            @views arho[:,1] .= rho_loc[:,1] .+ rho_core/2
            @views arho[:,2] .= rho_loc[:,2] .+ rho_core/2
            @views calc_epsxc_Vxc_VWN!( xc_calc, arho, e_rad, v_rad[ix] )
            #@views e_rad .= e_rad .* ( rho_rad[:,1] .+ rho_rad[:,2] .+ rho_core .* r2 )
            @views e_rad .= e_rad .* (arho[:,1] .+ arho[:,2]) .* r2
        end
    
        #println("sum e_rad = ", sum(e_rad))
        #println("sum v_rad[:,ix,1] = ", sum(v_rad[:,ix,1]))

        # Integrate to obtain the energy
        wx = pspotNL.paw.spheres[isp].ww[ix]
        ss = wx*PWDFT.integ_simpson( Nrmesh, e_rad, pspots[isp].rab )
        #println("integrated energy from ix = ", ss)
        energy += ss
  
    end

    #
    # Recompose the sph. harm. expansion
    lmax_loc = pspots[isp].lmax_rho + 1
    PAW_rad2lm!( ia, atoms, pspotNL, lmax_loc, v_rad, v_lm )

    return energy

end
