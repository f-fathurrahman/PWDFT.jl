function _driver_xc_PBE!(
    xc_calc::LibxcXCCalculator,
    Rhoe,
    gRhoe2,
    epsxc,
    V_xc,
    Vg_xc
)
    Npoints = size(Rhoe, 1)
    Nspin = 1

    eps_x = zeros(Float64,Npoints)
    eps_c = zeros(Float64,Npoints)

    V_x = zeros(Float64,Npoints)
    V_c = zeros(Float64,Npoints)
 
    Vg_x = zeros(Float64,Npoints)
    Vg_c = zeros(Float64,Npoints)
 
    ptr = Libxc_xc_func_alloc()
    # exchange part
    Libxc_xc_func_init(ptr, 101, Nspin)
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_gga_exc_vxc!(ptr, Npoints, Rhoe, gRhoe2, eps_x, V_x, Vg_x)
    Libxc_xc_func_end(ptr)

    #
    # correlation part
    Libxc_xc_func_init(ptr, 130, Nspin)
    Libxc_xc_func_set_dens_threshold(ptr, 1e-10)
    Libxc_xc_gga_exc_vxc!(ptr, Npoints, Rhoe, gRhoe2, eps_c, V_c, Vg_c)
    Libxc_xc_func_end(ptr)

    #
    Libxc_xc_func_free(ptr)

    # update the outputs:
    @views epsxc[:] .= eps_x[:] .+ eps_c[:]
    @views V_xc[:,1] .= V_x[:] .+ V_c[:]
    @views Vg_xc[:,1] = Vg_x[:] + Vg_c[:]  # for gradient correction

    # V_xc need to be corrected later using Vg_xc
    return
end


function _driver_xc_PBE!(
    xc_calc::XCCalculator,
    arho, grhoe2, sxc, v1xc, v2xc
)
    SMALL1 = 1e-6
    Nrmesh = size(arho, 1)
    for ir in 1:Nrmesh
        eex, vvx = XC_x_slater( arho[ir,1] )
        eec, vvc = XC_c_pw( arho[ir,1] )
        #
        sx, v1x, v2x = XC_x_pbe( arho[ir,1], grhoe2[ir] )
        sc, v1c, v2c = XC_c_pbe( arho[ir,1], grhoe2[ir] )
        #
        # NOTE: we are using Libxc convention
        # FIXME: not sure where to put the guard against small rho
        if arho[ir,1] > SMALL1
            sxc[ir] = (sx + sc)/arho[ir,1] + eex + eec # Libxc convention
        else
            sxc[ir] = eex + eec
        end
        v1xc[ir,1] = v1x + v1c + (vvx + vvc) # already the same as Libxc convention
        v2xc[ir,1] = (v2x + v2c)*0.5 # Libxc convention
    end
    return
end


# Use the density produced by sum_rad_rho to compute xc potential
# and energy, as xc functional is not diagonal on angular momentum
# numerical integration is performed.
function PAW_xc_potential_GGA!(
    AE::Bool, ia,
    atoms, pspots, pspotNL,
    xc_calc,
    rho_lm, v_lm
)

    isp = atoms.atm2species[ia] # species index
    Nrmesh = pspots[isp].Nr
    l2 = (pspots[isp].lmax_rho + 1)^2
    Nspin = size(rho_lm, 3)

    rho_rad = zeros(Float64, Nrmesh, Nspin)
    grad = zeros(Float64, Nrmesh, 3, Nspin) # gradient (r, ϕ, θ)
    grad2 = zeros(Float64, Nrmesh, Nspin) # square modulus of gradient

    nx = pspotNL.paw.spheres[isp].nx

    if AE
        rho_core = pspots[isp].paw_data.ae_rho_atc
    else
        rho_core = pspots[isp].rho_atc
    end

    r = pspots[isp].r
    r2 = r.^2

    arho = zeros(Float64, Nrmesh) # 2nd dim removed 
    grhoe2 = zeros(Float64, Nrmesh)


    # These arrays should depend on spin
    v1xc = zeros(Float64, Nrmesh, Nspin)
    v2xc = zeros(Float64, Nrmesh, Nspin)
    # These arrays don't depend on spin
    sxc = zeros(Float64, Nrmesh)

    # for energy, it will be summed over for all nx
    e_rad = zeros(Float64, Nrmesh)
    # for potential, will be processed later, depend on nx
    #gc_rad = zeros(Float64, Nrmesh, nx, Nspin)
    #h_rad = zeros(Float64, Nrmesh, 3, nx, Nspin)
    gc_rad = Vector{Matrix{Float64}}(undef, nx)
    h_rad = Vector{Array{Float64,3}}(undef, nx)
    for ix in 1:nx
        gc_rad[ix] = zeros(Float64, Nrmesh, Nspin)
        h_rad[ix] = zeros(Float64, Nrmesh, 3, Nspin)
    end

    energy = 0.0 # This should be accumulated for all ix
    spheres = pspotNL.paw.spheres

    #println("ia = ", ia, " AE = ", AE)
    #println("sum rho_core = ", sum(rho_core))
    #println("sum rho_lm = ", sum(rho_lm))

    for ix in 1:nx

        #println("Begin ix = ", ix)

        PAW_lm2rad!(ia, ix, atoms, pspots, pspotNL, rho_lm, rho_rad)
        #println("sum rho_rad = ", sum(rho_rad))
        
        PAW_gradient!(ia, ix, atoms, pspots, pspotNL,
            rho_lm, rho_rad, rho_core,
            grad2, grad
        )

        #@printf("%5d grad r    : %18.10e\n", ix, sum(grad[:,1,:]))
        #@printf("%5d grad phi  : %18.10e\n", ix, sum(grad[:,2,:]))
        #@printf("%5d grad theta: %18.10e\n", ix, sum(grad[:,3,:]))

        for ir in 1:Nrmesh
            arho[ir,1] = abs(rho_rad[ir,1]/r2[ir] + rho_core[ir])
            grhoe2[ir] = grad[ir,1]^2 + grad[ir,2]^2 + grad[ir,3]^2
        end

        _driver_xc_PBE!(xc_calc, arho, grhoe2, sxc, v1xc, v2xc)

        # radial stuffs
        # NOTE: We are using Libxc convention: e_rad is multiplied by arho and
        # h_rad will be multiplied by 2 (here or later in div_h, for calculating the
        # gradient correction to the potential)
        for ir in 1:Nrmesh
            e_rad[ir] = sxc[ir] * r2[ir] * arho[ir]
            gc_rad[ix][ir,1] = v1xc[ir,1]
            @views h_rad[ix][ir,1:3,1] = v2xc[ir,1]*grad[ir,1:3,1]*r2[ir]
            # Nspin == 1
        end
    
        # integrate to obtain the energy
        energy += integ_simpson(Nrmesh, e_rad, pspots[isp].rab)*spheres[isp].ww[ix]

        #ee = integ_simpson(Nrmesh, e_rad, pspots[isp].rab)*spheres[isp].ww[ix]
        #println("energy for current ix = ", ix, " ", ee)

    end

    #println("sum gc_rad = ", sum(gc_rad))
    #println("sum h_rad = ", sum(h_rad))

    lmax_loc = pspots[isp].lmax_rho + 1
    gc_lm = zeros(Float64, Nrmesh, l2, Nspin)
    # convert the first part of the GC correction back to spherical harmonics
    PAW_rad2lm!( ia, atoms, pspotNL, lmax_loc, gc_rad, gc_lm)
    #println("sum gc_lm = ", sum(gc_lm))

    # trick to get faster convergence w.r.t to θ
    for ix in 1:nx
        @views h_rad[ix][:,3,1] = h_rad[ix][:,3,1] / spheres[isp].sin_th[ix]
        # Nspin = 1
    end
    # 
    #println("sum h_rad after divided by sin_th = ", sum(h_rad))

    # We need the gradient of H to calculate the last part of the exchange
    # and correlation potential. First we have to convert H to its Y_lm expansion
    #PAW_rad2lm!( i, h_rad, h_lm, i%l+rad(i%t)%ladd, nspin_gga )
    lmax_loc_add = lmax_loc + spheres[isp].ladd
    h_lm = zeros(Float64, Nrmesh, 3, lmax_loc_add^2, Nspin)

    #println("lmax_loc = ", lmax_loc)
    #println("ladd = ", spheres[isp].ladd)
    #println("lmax_loc_add = ", lmax_loc_add)
    
    #@views PAW_rad2lm!( ia, atoms, pspotNL, lmax_loc_add, h_rad[:,1,:,:], h_lm[:,1,:,:])
    #@views PAW_rad2lm!( ia, atoms, pspotNL, lmax_loc_add, h_rad[:,2,:,:], h_lm[:,2,:,:])
    #@views PAW_rad2lm!( ia, atoms, pspotNL, lmax_loc_add, h_rad[:,3,:,:], h_lm[:,3,:,:])
    PAW_rad2lm3!(ia, atoms, pspotNL, lmax_loc_add, h_rad, h_lm)
    #println("sum h_lm = ", sum(h_lm))

    # Calculate the divergence: ∇ ⋅ h_rad = div_h
    div_h = zeros(Float64, Nrmesh, lmax_loc^2, Nspin)
    PAW_divergence!(
        ia, atoms, pspots, pspotNL,
        h_lm, div_h, lmax_loc_add, lmax_loc
    )
    #println("sum div_h = ", sum(div_h))

    #println("sum abs v_lm before in PAW_xc_potential_GGA: ", sum(abs.(v_lm)))

    # Finally sum it back into v_xc
    # Factor 2 of div_h because we are using Libxc convention
    for ispin in 1:Nspin
        for lm in 1:l2
            @views v_lm[1:Nrmesh,lm,ispin] .= gc_lm[1:Nrmesh,lm,ispin] .- 2*div_h[1:Nrmesh,lm,ispin]
        end
    end

    #println("energy in PAW_xc_potential_GGA: ", energy)
    
    #println("size v_lm = ", size(v_lm))
    #println("sum abs v_lm after in PAW_xc_potential_GGA: ", sum(abs.(v_lm)))

    # Energy is scalar, it is returned
    return energy

end