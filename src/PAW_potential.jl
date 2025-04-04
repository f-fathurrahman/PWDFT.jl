function PAW_potential!(Ham::Hamiltonian)
    return PAW_potential!(
        Ham.atoms, Ham.pspots, Ham.pspotNL, Ham.xc_calc,
        Ham.pspotNL.paw.becsum, # XXX FIXME should be rho_bec (symmetrized and only for PAW)
        Ham.pspotNL.paw.ddd_paw,
        Ham.pspotNL.paw.E_paw_cmp
    )
end


function PAW_potential!(
    atoms::Atoms,
    pspots::Vector{PsPot_UPF},
    pspotNL::PsPotNL_UPF,
    xc_calc,
    becsum, ddd_paw, e_cmp
)

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
    # This is needed to minimize arrays allocation; Nrmesh is generally
    # different for different species
    # XXX: Alternative: use Nrmesh_max

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

        #println("\nLoop PAW_potential: ia = ", ia)

        isp = atm2species[ia]
        # Skip this atoms if it is not using PAW
        if !pspots[isp].is_paw
            continue
        end

        Nrmesh = pspots[isp].Nr
        l2 = (pspots[isp].lmax_rho + 1)^2
        kkbeta = pspots[isp].kkbeta

        # We might need to call GC manually to free these arrays (if
        # they are allocated inside the loop).
        # For the moment we allocate per-species arrays
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

            # Compute rho_lm from becsum
            PAW_rho_lm!(AE, ia, atoms, pspots, pspotNL, becsum, rho_lm)
            #
            #println("\nsum becsum total = ", sum(becsum))
            #println("sum becsum spin up = ", sum(becsum[:,:,1]))
            #(Nspin == 2) && println("sum becsum spin dn = ", sum(becsum[:,:,2]))
            #println("AE = $(AE) sum rho_lm total = $(sum(rho_lm))")
            #println("AE = $(AE) sum rho_lm up = $(sum(rho_lm[:,:,1]))")
            #(Nspin == 2) && println("AE = $(AE) sum rho_lm dn = $(sum(rho_lm[:,:,2]))")

            # Hartree term
            #println("\nCalling PAW_h_potential")
            @views energy = PAW_h_potential!( ia, atoms, pspots, rho_lm, v_lm[:,:,1] )
            #println("After calling PAW_h_potential: energy = ", energy)
            energy_tot += sgn*energy
            e_cmp[ia,1,i_what] = sgn*energy # Hartree, all-electron
            for ispin in 1:Nspin # ... v_H has to be copied to all spin components
                @views savedv_lm[:,:,ispin] = v_lm[:,:,1]
            end

            # XC term
            #println("\nCalling PAW_xc_potential")
            if xc_calc.family == :GGA
                energy = PAW_xc_potential_GGA!( AE, ia, atoms, pspots, pspotNL, xc_calc, rho_lm, v_lm )
            else
                # Default is LDA
                energy = PAW_xc_potential_LDA!( AE, ia, atoms, pspots, pspotNL, xc_calc, rho_lm, v_lm )
            end
            # FIXME: metaGGA is not yet supported
            #println("After calling PAW_xc_potential: energy = ", energy)
            #println("After calling PAW_xc_potential: sum v_lm up = ", sum(v_lm[:,:,1]))
            #(Nspin == 2) && println("After calling PAW_xc_potential: sum v_lm dn = ", sum(v_lm[:,:,2]))

            energy_tot += sgn*energy
            e_cmp[ia,2,i_what] = sgn*energy # XC, all-electron
            # PAW_xc_func! is alias to either PAW_xc_potential! or PAW_xc_potential_GGA!
            # switch between LDA and GGA are made above

            @views savedv_lm[:,:,1:Nspin] .+= v_lm[:,:,1:Nspin]
            
            #println("\nBefore calc ddd_paw: sum savedv_lm up = ", sum(savedv_lm[:,:,1]))
            #(Nspin == 2) && println("Before calc ddd_paw: sum savedv_lm dn = ", sum(savedv_lm[:,:,2]))

            #
            # Calculate ddd_paw
            #
            for ispin in 1:Nspin
                nmb = 0
                # loop on all pfunc for this kind of pseudo
                for nb in 1:nh[isp], mb in nb:nh[isp]
                    nmb = nmb + 1
                    #
                    # compute the density from a single pfunc
                    becfake[nmb,ia,ispin] = 1.0
                    PAW_rho_lm!(AE, ia, atoms, pspots, pspotNL, becfake, rho_lm)
                    #println("nmb = $(nmb) sum rho_lm spin up = ", sum(rho_lm[:,:,1]))
                    #println("nmb = $(nmb) sum rho_lm spin dn = ", sum(rho_lm[:,:,2]))
                    #
                    # Now I multiply the rho_lm and the potential, I can use
                    # rho_lm itself as workspace
                    #
                    for lm in 1:l2
                        for ir in 1:Nrmesh
                            rho_lm[ir,lm,ispin] *= savedv_lm[ir,lm,ispin]
                        end
                        # Integrate (XXX: why using kkbeta?)
                        #println("Nrmesh = ", Nrmesh, " kkbeta = ", kkbeta)
                        @views res = PWDFT.integ_simpson( kkbeta, rho_lm[:,lm,ispin], pspots[isp].rab )
                        ddd_paw[nmb,ia,ispin] += sgn*res
                    end
                    # restore becfake to zero
                    becfake[nmb,ia,ispin] = 0.0
                end
            end # Nspin
        end # AE, PS

        #println("\nia=$(ia) sum ddd_paw up = $(sum(ddd_paw[:,ia,1]))")
        #(Nspin == 2) && println("ia=$(ia) sum ddd_paw dn = $(sum(ddd_paw[:,ia,2]))")

    end # loop over all atoms

    return energy_tot
end
