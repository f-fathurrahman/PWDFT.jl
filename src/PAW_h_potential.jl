function PAW_h_potential!(
    ia,
    atoms::Atoms, pspots,
    rho_lm,
    v_lm
)
    #=
    !---------------------------------------------------------------------------------
    !! Computes H  potential from rho, used by PAW_h_energy and PAW_ddot.
    !
    USE constants,              ONLY : fpi, e2
    USE radial_grids,           ONLY : hartree
    USE uspp_param,             ONLY : upf
    USE noncollin_module,       ONLY : nspin_lsda
    USE ions_base,              ONLY : ityp
    USE lsda_mod,               ONLY : nspin
    USE atom,                   ONLY : g => rgrid
    !
    TYPE(paw_info), INTENT(IN) :: i !! atom's minimal info
    REAL(DP), INTENT(IN)  :: rho_lm(i%m,i%l**2,nspin) !! charge density as lm components already summed on spin
    REAL(DP), INTENT(OUT) :: v_lm(i%m,i%l**2) !! potential as lm components
    !
    ! ... local variables
    !
    REAL(DP) :: aux(i%m) ! workspace
    REAL(DP) :: pref     ! workspace
    !
    INTEGER  :: lm,l     ! counter on composite angmom lm = l**2 +m
    INTEGER  :: k        ! counter on radial grid (only for energy) 
    REAL(DP) :: e        ! workspace
    =#
    
    #
    # this loop computes the hartree potential using the following formula:
    #               l is the first argument in hartree subroutine
    #               r1 = min(r,r'); r2 = MAX(r,r')
    # V_h(r) = \sum{lm} Y_{lm}(\hat{r})/(2l+1) \int dr' 4\pi r'^2 \rho^{lm}(r') (r1^l/r2^{l+1})
    #     done here --> ^^^^^^^^^^^^^^^^^^^^^           ^^^^^^^^^^^^^^^^^^^^^^ <-- input to the hartree subroutine
    #                 output from the h.s. --> ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    atm2species = atoms.atm2species

    isp = atm2species[ia]
    l2 = (pspots[isp].lmax_rho + 1)^2
    r = pspots[isp].r
    dx = pspots[isp].dx

    Nrmesh = pspots[isp].Nr
    aux = zeros(Float64, Nrmesh)
    nspin_lsda = size(rho_lm, 3)

    fill!(v_lm, 0.0)

    for lm in 1:l2
        l = floor(Int64, sqrt(lm - 1)) # l has to start from *zero*
        pref = 4*pi/(2*l + 1)
        for k in 1:Nrmesh
            aux[k] = pref * sum(rho_lm[k,lm,1:nspin_lsda])
        end
        @views radial_hartree!( l, 2*l + 2, r, dx, aux, v_lm[:,lm] )
    end
    #
    # compute energy (we compute here by default):
    # E_h = \sum_lm \int v_lm(r) (rho_lm(r) r^2) dr
    energy = 0.0
    for lm in 1:l2
        for k in 1:Nrmesh
            aux[k] = v_lm[k,lm] * sum(rho_lm[k,lm,1:nspin_lsda])
        end
        ene = PWDFT.integ_simpson( Nrmesh, aux, pspots[isp].rab )
        # Sum all the energies in PAW_ddot
        energy += ene
    end
    # fix double counting
    energy *= 0.5

    return energy
end