function PAW_gradient!(
    ia, ix,
    atoms, pspots, pspotNL,
    rho_lm, rho_rad, rho_core, grho_rad2, grho_rad
)
    # Build gradient of radial charge distribution from its spherical harmonics expansion

#=
    INTEGER, INTENT(IN)  :: ix
    !! line of the dylm2 matrix to use actually it is  one of the nx spherical
    !! integration directions
    
    - rho_lm(i%m,i%l**2,nspin_gga) !! Y_lm expansion of rho
    - rho_rad(i%m,nspin_gga) !! radial density along direction ix
    - rho_core(i%m) !! core density
    REAL(DP), INTENT(OUT) :: grho_rad2(i%m,nspin_gga) !! |grad(rho)|^2 on rad. grid
    REAL(DP), OPTIONAL,INTENT(OUT):: grho_rad(i%m,3,nspin_gga) !! vector gradient (only for gcxc)
=#

    atm2species = atoms.atm2species
    isp = atm2species[ia]
    Nrmesh = pspots[isp].Nr

    Nspin = size(rho_lm, 3)

    # r, theta and phi components ---^
    #
    #REAL(DP) :: aux(i%m), aux2(i%m), fact  ! workspace
    aux = zeros(Float64, Nrmesh)
    aux2 = zeros(Float64, Nrmesh)

    # 1. build real charge density = rho/r**2 + rho_core
    # 2. compute the partial derivative of rho_rad
    fact = 1.0/Nspin
    fill!(grho_rad2, 0.0)
    for ispin in 1:Nspin
        # build real charge density
        @views aux[1:Nrmesh] .= rho_rad[1:Nrmesh,ispin] ./ r2[1:Nrmesh] .+ rho_core[1:Nrmesh]*fact
        radial_gradient_AE( r, aux, aux2 ) # aux2 is the output
        # compute the square
        @views grho_rad2[:,ispin] .= aux2[:].^2
        # store in vector gradient, if present:
        @views grho_rad[:,1,is] .= aux2[:]
    end
 
    for ispin in 1:Nspin
        fill!(aux, 0.0)
        fill!(aux2, 0.0)
        # Spherical (lm=1) component (that would also include core correction) can be omitted
        # as its contribution to non-radial derivative is zero
        for lm in 2:lm_max
            # 5. [ \sum_{lm} rho(r) (dY_{lm}/dphi /cos(theta))  ]**2
            aux[1:Nrmesh] .= aux[1:Nrmesh] .+ dylmp[ix,lm] * rho_lm[1:Nrmesh,lm,ispin]
            # 6. [ \sum_{lm} rho(r) (dY_{lm}/dtheta)  ]**2
            aux2[1:Nrmesh] .= aux2[1:Nrmesh] .+ dylmt[ix,lm] * rho_lm[1:Nrmesh,lm,ispin]
        end
        # Square and sum up these 2 components, the (1/r**2)**3 factor come from:
        #  a. 1/r**2 from the derivative in spherical coordinates
        #  b. (1/r**2)**2 from rho_lm being multiplied by r**2 
        #     (as the derivative is orthogonal to r you can multiply after deriving)
        @views grho_rad2[1:Nrmesh,ispin] .+= ( aux[1:Nrmesh]^2 .+ aux2[1:Nrmesh]^2) / r[:].^6
        # Store vector components:
        @views grho_rad[1:Nrmesh,2,ispin] .= aux[1:Nrmesh]  ./ r[:].^3 # phi 
        @views grho_rad[1:Nrmesh,3,ispin] .= aux2[1:Nrmesh] ./ r[:].^3  # theta
    end

    return
end