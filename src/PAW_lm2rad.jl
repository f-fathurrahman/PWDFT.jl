# Build radial charge distribution from its spherical harmonics expansion.
function PAW_lm2rad!(
    ia, ix,
    atoms, pspots, pspotNL,
    F_lm, F_rad
)
    #INTEGER :: ix
    # !! line of the ylm matrix to use
    # !! actually it is one of the nx directions
  
    #INTEGER, INTENT(IN) :: nspin  #!! number of spin components
    #REAL(DP), INTENT(IN) :: F_lm(i%m,i%l**2,nspin) # Y_lm expansion of rho
    #REAL(DP), INTENT(OUT) :: F_rad(i%m,nspin)

    Nspin = size(F_rad,2)
    isp = atoms.atm2species[ia]
    l2 = (pspots[isp].lmax_rho + 1)^2
    sphere = pspotNL.paw.spheres[isp]

    fill!(F_rad, 0.0)
    for ispin in 1:Nspin, lm in 1:l2
        F_rad[:,ispin] .+= sphere.ylm[ix,lm]*F_lm[:,lm,ispin]
    end
    return
end