function PAW_potential!(
    atoms, pspots, pspotNL, xc_calc,
    becsum, ddd_paw, e_cmp
)

    # ddd_paw: descreening coefficients (AE - PS)
    # REAL(DP), INTENT(IN) :: becsum(nhm*(nhm+1)/2,nat,nspin) !! cross band occupations
    # REAL(DP), INTENT(OUT) :: d(nhm*(nhm+1)/2,nat,nspin) !! descreening coefficients (AE - PS)
    # REAL(DP), INTENT(OUT), OPTIONAL :: e_cmp(nat,2,2) !! components of the energy

    becfake = zeros(Float64, size(becsum))
    fill!(ddd_paw, 0.0)
    fill!(e_cmp, 0.0)
    energy_tot = 0.0

    Nspin = size(ddd_paw, 3)

    Nspecies = atoms.Nspecies
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    nh = pspotNL.nh

    #Nrmesh_max = 0
    #for isp in 1:atoms.Nspecies
    #    if Nrmesh_max < pspots[isp].Nr
    #        Nrmesh_max = Nr
    #    end
    #end
    #println("Nrmesh_max = ", Nrmesh_max)
    #v_lm = zeros(Float64, Nrmesh_max, l2, Nspin)
    #savedv_lm = zeros(Float64, Nrmesh_max, l2, Nspin)
    #rho_lm = zeros(Float64, Nrmesh_max, l2, Nspin)

    # Preallocate work arrays (for each species)
    v_lm_s = Vector{Array{Float64,3}}(undef,atoms.Nspecies)
    savedv_lm_s = Vector{Array{Float64,3}}(undef,atoms.Nspecies)
    rho_lm_s = Vector{Array{Float64,3}}(undef,atoms.Nspecies)
    for isp in 1:Nspecies
        Nrmesh = pspots[isp].Nr
        l2 = (pspots[isp].lmax_rho + 1)^2
        v_lm_s[isp] = zeros(Float64, Nrmesh, l2, Nspin)
        savedv_lm_s[isp] = zeros(Float64, Nrmesh, l2, Nspin)
        rho_lm_s[isp] = zeros(Float64, Nrmesh, l2, Nspin)
    end


    # Begin loop over all atoms
    for ia in 1:Natoms

        isp = atm2species[ia]
        # Skip this atoms if it is not using PAW
        if !pspots[isp].is_paw
            continue
        end

        Nrmesh = pspots[isp].Nr
        l2 = (pspots[isp].lmax_rho + 1)^2
        kkbeta = pspots[isp].kkbeta

        # We might need to call GC manually to free these arrays
        v_lm = v_lm_s[isp]
        savedv_lm = savedv_lm_s[isp]
        rho_lm = rho_lm_s[isp]

        # These arrays will be set to zeros in the corresponding caller functions
        # savedv_lm will be reset here
        fill!(savedv_lm, 0.0)

        # All-electron contribution
        for AE in [true, false]
            if AE
                i_what = 1
                sgn = 1
            else
                i_what = 2
                sgn = -1
            end
            # sgn: sign for energy summation

            PAW_rho_lm!(AE, ia, atoms, pspots, pspotNL, becsum, rho_lm)

            # Hartree term
            @views energy = PAW_h_potential!( ia, atoms, pspots, rho_lm, v_lm[:,:,1] )
            energy_tot += sgn*energy
            e_cmp[ia,1,i_what] = sgn*energy # Hartree, all-electron
            for ispin in 1:Nspin # ... v_H has to be copied to all spin components
                @views savedv_lm[:,:,ispin] = v_lm[:,:,1]
            end

            # XC term
            @views energy = PAW_xc_potential!( AE, ia, atoms, pspots, pspotNL, xc_calc, rho_lm, v_lm )
            energy_tot += sgn*energy
            e_cmp[ia,2,i_what] = sgn*energy # XC, all-electron

            savedv_lm[:,:,:] .+= v_lm[:,:,:]

            for ispin in 1:Nspin
                nmb = 0
                # loop on all pfunc for this kind of pseudo
                for nb in 1:nh[isp], mb in nb:nh[isp]
                    nmb = nmb + 1 # nmb = 1, nh*(nh+1)/2
                    #
                    # compute the density from a single pfunc
                    becfake[nmb,ia,ispin] = 1.0
                    PAW_rho_lm!(AE, ia, atoms, pspots, pspotNL, becfake, rho_lm)
                    #
                    # Now I multiply the rho_lm and the potential, I can use
                    # rho_lm itself as workspace
                    #
                    for lm in 1:l2
                        for ir in 1:Nrmesh
                            rho_lm[ir,lm,ispin] *= savedv_lm[ir,lm,ispin]
                        end
                        # Integrate (XXX: why using kkbeta?)
                        res = PWDFT.integ_simpson( kkbeta, rho_lm[:,lm,ispin], pspots[isp].rab )
                        ddd_paw[nmb,ia,ispin] += sgn*res
                    end
                    # restore becfake to zero
                    becfake[nmb,ia,ispin] = 0.0
                end
            end # Nspin

        end # AE, PS
    end # loop over all atoms

    return energy_tot
end