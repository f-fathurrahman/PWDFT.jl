function PAW_gcxc_potential( # i, rho_lm, rho_core, v_lm, energy )
    ia,
    atoms::Atoms, pspots, pspotNL,
    xc_calc,
    rho_lm, rho_core, v_lm, energy
)
    # Add gradient correction to v_xc, code mostly adapted from ../atomic/vxcgc.f90
    # in order to support non-spherical charges (as Y_lm expansion).  
    # Note that the first derivative in vxcgc becomes a gradient, while the second is
    # a divergence.  
    # We also have to temporarily store some additional Y_lm components in order not
    # to loose precision during the calculation, even if only the ones up to 
    # lmax_rho (the maximum in the density of charge) matter when computing \int v*rho.

#=
    !
    TYPE(paw_info), INTENT(IN) :: i
    !! atom's minimal info
    REAL(DP), INTENT(IN) :: rho_lm(i%m,i%l**2,nspin)
    !! charge density as lm components
    REAL(DP), INTENT(IN) :: rho_core(i%m)
    !! core charge, radial and spherical
    REAL(DP), INTENT(INOUT) :: v_lm(i%m,i%l**2,nspin)
    !! potential to be updated
    REAL(DP), OPTIONAL, INTENT(INOUT) :: energy
    !! if present, add GC to energy
    !
    ! ... local variables
    !
    REAL(DP), PARAMETER :: epsr = 1.e-6_DP, epsg = 1.e-10_DP
    ! (as in PW/src/gradcorr.f90)
    !
    REAL(DP), ALLOCATABLE :: rho_rad(:,:) ! charge density sampled
    REAL(DP), ALLOCATABLE :: grad(:,:,:)  ! gradient
    REAL(DP), ALLOCATABLE :: gradx(:,:,:) ! gradient (swapped indexes)
    REAL(DP), ALLOCATABLE :: grad2(:,:)   ! square modulus of gradient
                                          ! (first of charge, than of hamiltonian)
    REAL(DP), ALLOCATABLE :: gc_rad(:,:,:)    ! GC correction to V (radial samples)
    REAL(DP), ALLOCATABLE :: gc_lm(:,:,:)     ! GC correction to V (Y_lm expansion)
    REAL(DP), ALLOCATABLE :: h_rad(:,:,:,:)   ! hamiltonian (vector field)
    REAL(DP), ALLOCATABLE :: h_lm(:,:,:,:)    ! hamiltonian (vector field)
                                    ! ^^^^^^^^^^^^^^^^^^ expanded to higher lm than rho !
    REAL(DP), ALLOCATABLE :: div_h(:,:,:)  ! div(hamiltonian)
    !
    REAL(DP), ALLOCATABLE :: rhoout_lm(:,:,:) ! charge density as lm components
    REAL(DP), ALLOCATABLE :: vout_lm(:,:,:)   ! potential as lm components
    REAL(DP), ALLOCATABLE :: segni_rad(:,:)   ! sign of the magnetization
    !
    REAL(DP), ALLOCATABLE :: arho(:,:), grad2_v(:)
    REAL(DP), ALLOCATABLE :: r_vec(:,:)
    !
    REAL(DP), DIMENSION(i%m,nspin_gga) :: v1x, v2x, v1c, v2c  !workspace
    REAL(DP), DIMENSION(i%m) :: sx, sc
    REAL(DP), ALLOCATABLE :: v2cud(:)
    !
    REAL(DP) :: vnull
    !
    REAL(DP), ALLOCATABLE :: e_rad(:)      ! aux, used to store energy
    REAL(DP) :: e, e_gcxc                  ! aux, used to integrate energy
    !
    INTEGER  :: k, ix, is, lm              ! counters on spin and mesh
    REAL(DP) :: sgn                        ! workspace
    REAL(DP) :: co2                        ! workspace
    !
    INTEGER :: mytid, ntids

    REAL(DP),ALLOCATABLE :: egcxc_of_tid(:)

=#
  
    atm2species = atoms.atm2species
    isp = atm2species[ia]
    Nrmesh = pspots[isp].Nr

    Nspin = size(rho_lm, 3) # nspin_gga

    @assert Nspin == 1 # current limitation

    spheres = pspotNL.paw.spheres
    nx = spheres[isp].nx

    e_gcxc = 0.0

    lm_max = (spheres[isp].lmax + 1)^2
    gc_rad = zeros(Float64, Nrmesh, nx, Nspin) # GC correction to V (radial samples)
    gc_lm = zeros(Float64, Nrmesh, lm_max, Nspin) # GC correction to V (Y_lm expansion)
    h_rad_x = zeros(Float64, Nrmesh, nx, Nspin) # hamiltonian (vector field)
    h_rad_y = zeros(Float64, Nrmesh, nx, Nspin) # hamiltonian (vector field)
    h_rad_z = zeros(Float64, Nrmesh, nx, Nspin) # hamiltonian (vector field)

    lm_max_add = (spheres[isp].lmax + 1 + spheres[isp].ladd)^2
    h_lm = zeros(Float64, Nr, lm_max_add, Nspin)
    # expanded to higher lm than rho

    div_h = zeros(Float64, Nrmesh, lm_max, Nspin)

    # charge density as lm components
    rhoout_lm = zeros(Float64, Nrmesh, lm_max, Nspin)

    # potential as lm components
    vout_lm = zeros(Float64, Nrmesh, lm_max, Nspin)

    # charge density as lm components
    segni_rad = zeros(Float64, Nrmesh, nx)      
    
    if Nspin == 2
        # XXX: Non-collinear case is removed
        rhoout_lm[:] .= rho_lm[:]
    end

    rho_rad = zeros(Float64, Nrmesh, Nspin) # charge density sampled
    grad = zeros(Float64, Nrmesh, 3, Nspin) # gradient
    
    # square modulus of gradient
    # (first of charge, then of hamiltonian)
    grad2 = zeros(Float64, Nrmesh, Nspin)

    e_rad = zeros(Float64, Nrmesh)
    

    # Nspin == 1 case

    arho = zeros(Float64, Nrmesh, 1) # XXX why need 2nd dim ?
    grad2_v = zeros(Float64, Nrmesh)
    gradx = zeros(Float64, 3, Nrmesh, 1) # XXX why need 2nd dim ?
    # May be because PAW_lm2rad need 2d array
    
    for ix in 1:nx
        #
        #  WARNING: the next 2 calls are duplicated for spin==2
        PAW_lm2rad!( ia, ix, rho_lm, rho_rad, nspin_mag )
        PAW_gradient!( ia, ix, rho_lm, rho_rad, rho_core, grad2, grad )
        
        for k in 1:Nrmesh
            arho[k,1] = rho_rad[k,1]/r2[k] + rho_core[k]
            arho[k,1] = abs(arho[k,1])
            @views gradx[1:3,k,1] .= grad[k,1:3,1]
        end
        
        #CALL xc_gcx( i%m, 1, arho, gradx, sx, sc, v1x, v2x, v1c, v2c )
        # CALL GGA subroutine here
        
        # radial stuffs
        for k in 1:Nrmesh
            e_rad[k] = (sx[k] + sc[k]) * r2[k]
            gc_rad[k,ix,1]  = ( v1x[k,1] + v1c[k,1] )  # * g(i%t)%rm2(k)
            h_rad[k,1:3,ix,1] = ( v2x[k,1] + v2c[k,1] )*grad[k,1:3,1]*r2[k]
        end
        
        # integrate to obtain the energy
        energy += integ_simpson(Nrmesh, e_rad, rab)*spheres[isp].ww[ix]

    end
        
    # DEALLOCATE( arho, grad2_v ) 
    # DEALLOCATE( gradx )

#=

    ELSEIF ( nspin_mag == 2 .OR. nspin_mag == 4 ) THEN
        !
        ALLOCATE( gradx(3,i%m,2) )
        ALLOCATE( r_vec(i%m,2) )
        ALLOCATE( v2cud(i%m) )
        !
        !   this is the \sigma-GGA case
        !
        DO ix = ix_s, ix_e
           !
           CALL PAW_lm2rad( i, ix, rhoout_lm, rho_rad, nspin_gga )
           CALL PAW_gradient( i, ix, rhoout_lm, rho_rad, rho_core,grad2, grad )
           !
           DO k = 1, i%m
               !
               ! Prepare the necessary quantities
               ! rho_core is considered half spin up and half spin down:
               co2 = rho_core(k)/2
               ! than I build the real charge dividing by r**2
               r_vec(k,1) = rho_rad(k,1)*g(i%t)%rm2(k) + co2
               r_vec(k,2) = rho_rad(k,2)*g(i%t)%rm2(k) + co2
               !
               !
               gradx(:,k,1) = grad(k,:,1)
               gradx(:,k,2) = grad(k,:,2)
           ENDDO
           !
           CALL xc_gcx( i%m, 2, r_vec, gradx, sx, sc, v1x, v2x, v1c, v2c, v2cud )
           !
           DO k = 1, i%m
              !
              IF ( PRESENT(energy) ) e_rad(k) = e2*(sx(k)+sc(k))*g(i%t)%r2(k)
              !
              ! first term of the gradient correction : D(rho*Exc)/D(rho)
              gc_rad(k,ix,1)  = (v1x(k,1)+v1c(k,1)) !*g(i%t)%rm2(k)
              gc_rad(k,ix,2)  = (v1x(k,2)+v1c(k,2)) !*g(i%t)%rm2(k)
              !
              ! h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
              ! h_rad(k,:,ix,1) =( (v2xup_vec(k)+v2c)*grad(k,:,1)+v2c*grad(k,:,2) )*g(i%t)%r2(k)
              ! h_rad(k,:,ix,2) =( (v2xdw_vec(k)+v2c)*grad(k,:,2)+v2c*grad(k,:,1) )*g(i%t)%r2(k)
              h_rad(k,:,ix,1) =( (v2x(k,1)+v2c(k,1))*grad(k,:,1) + &
                                  v2cud(k)*grad(k,:,2) )*g(i%t)%r2(k)
              h_rad(k,:,ix,2) =( (v2x(k,2)+v2c(k,2))*grad(k,:,2) + &
                                  v2cud(k)*grad(k,:,1) )*g(i%t)%r2(k)
              !
           ENDDO
           !
           ! integrate energy (if required)
           ! NOTE: this integration is duplicated for every spin, FIXME!
           IF (PRESENT(energy)) THEN
               CALL simpson( i%m, e_rad, g(i%t)%rab, e )
               egcxc_of_tid(mytid) = egcxc_of_tid(mytid) + e * rad(i%t)%ww(ix)
           ENDIF
           !
        ENDDO ! ix
        !
        DEALLOCATE( gradx )
        DEALLOCATE( r_vec )
        DEALLOCATE( v2cud )
        !
    ELSE spin
    !
    ENDIF spin
    !
=#


    # convert the first part of the GC correction back to spherical harmonics
    PAW_rad2lm!( ia, gc_rad, gc_lm, i%l, nspin_gga )
    #
    # Note that the expansion into spherical harmonics of the derivative 
    # with respect to theta of the spherical harmonics, is very slow to
    # converge and would require a huge angular momentum ladd.
    # This derivative divided by sin_th is much faster to converge, so
    # we divide here before calculating h_lm and keep into account for
    # this factor sin_th in the expression of the divergence.
    #
    # ADC 30/04/2009.
    # 
    for ix in 1:nx
        h_rad[1:Nrmesh,3,ix,1:Nspin] = h_rad[1:Nr,3,ix,1:Nspin] / spheres[isp].sin_th[ix]
    end

    # We need the gradient of H to calculate the last part of the exchange
    # and correlation potential. First we have to convert H to its Y_lm expansion
    PAW_rad2lm3( ia, h_rad, h_lm, lm_max_add, nspin_gga )
    #
    # Compute div(H)
    PAW_divergence!( ia, h_lm, div_h, lm_max_add, i%l )
    #                       input max lm --^  output max lm-^
    

    # Finally sum it back into v_xc
    for ispin in 1:Nspin
        for lm in 1:lm_max
            @views vout_lm[1:Nrmesh,lm,ispin] .+= ( gc_lm[1:Nrmesh,lm,ispin] .- div_h[1:Nrmesh,lm,ispin] )
        end
    end

    # Noncollinear stuffs (with Nspin=4) are skipped

    v_lm[:,:,1:Nspin] .+= vout_lm[:,:,1:Nspin]

    return

end