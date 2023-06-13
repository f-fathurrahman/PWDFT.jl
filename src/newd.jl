function calc_integ_QVeff!( Ham )

    Nspin = Ham.electrons.Nspin
    Ng = Ham.pw.gvec.Ng
    lmaxkb = Ham.pspotNL.lmaxkb
    G2 = Ham.pw.gvec.G2
    G = Ham.pw.gvec.G
    Npoints = prod(Ham.pw.Ns)
    idx_g2r = Ham.pw.gvec.idx_g2r
    Veff = Ham.potentials.Total

    lmaxq = 2*lmaxkb + 1 # using 1-indexing
    ylmk0 = zeros(Float64, Ng, lmaxq*lmaxq)
    _lmax = lmaxq - 1 # or 2*lmaxkb
    # Ylm_real_qe accept l value starting from 0
    Ylm_real_qe!(_lmax, G, ylmk0)

    VeffG = zeros(ComplexF64, Ng, Nspin)

    # Fourier transform of the total effective potential
    ctmp = zeros(ComplexF64, Npoints)
    for ispin in 1:Nspin
        for ip in 1:Npoints
            ctmp[ip] = Veff[ip,ispin] # Veff already contains Ps_loc
        end
        #println("integ VQeff: sum Veff = ", sum(ctmp))
        R_to_G!(Ham.pw, ctmp)
        #ctmp /= sqrt(Npoints) # XXX: scale
        for ig in 1:Ng
            ip = idx_g2r[ig]
            VeffG[ig,ispin] = ctmp[ip]/Npoints
            # (1/Npoints) is needed for QE's convention
        end
        #println("sum VeffG = ", sum(VeffG))
    end


    Nspecies = Ham.atoms.Nspecies
    Natoms = Ham.atoms.Natoms
    atm2species = Ham.atoms.atm2species
    atpos = Ham.atoms.positions
    nh = Ham.pspotNL.nh
    Deeq = Ham.pspotNL.Deeq
    CellVolume = Ham.pw.CellVolume

    # Deeq will be recalculated here
    # No need to reset it to zeros.

    for isp in 1:Nspecies

        psp = Ham.pspots[isp]
        need_augmentation = psp.is_ultrasoft || psp.is_paw
        if !need_augmentation
            continue
        end

        # nij = max number of (ih,jh) pairs per atom type nt
        nij = round(Int64, nh[isp] * (nh[isp] + 1)/2 )

        Qgm = zeros(ComplexF64, Ng, nij)
        #println("nij = ", nij)
        #
        # Compute and store Q(G) for this atomic species 
        # (without structure factor)
        ijh = 0
        for ih in 1:nh[isp], jh in ih:nh[isp]
            ijh = ijh + 1
            #qvan2!( ngm, ih, jh, nt, qmod, qgm(1,ijh), ylmk0 )
            @views qvan2!( Ham.pspotNL, ih, jh, isp, G2, ylmk0, Qgm[:,ijh] )
            #sQ = sum(Qgm[:,ijh])
            #@printf("ijh = %4d sum Qgm = (%18.10f,%18.10f)\n", ijh, real(sQ), imag(sQ))
        end
        #
        # count max number of atoms of type isp
        nab = sum(atm2species .== isp)
        #
        aux = zeros(ComplexF64, Ng, nab)
        #
        # Compute and store V(G) times the structure factor e^(-iG*tau)
        for ispin in 1:Nspin
            nb = 0
            for ia in 1:Natoms
                if atm2species[ia] != isp
                    continue
                end
                nb = nb + 1
                for ig in 1:Ng
                    GX = atpos[1,ia]*G[1,ig] + atpos[2,ia]*G[2,ig] + atpos[3,ia]*G[3,ig]
                    Sf = cos(GX) + im*sin(GX) # conjugate
                    aux[ig,nb] = VeffG[ig,ispin] * Sf
                end
                #println("nb = ", nb, " sum aux VeffG*Sf = ", sum(aux[:,nb]))
            end
            #
            # here we compute the integral Q*V for all atoms of this kind
            deeaux = real(Qgm' * aux)
            #println("size deeaux = ", size(deeaux))
            #CALL DGEMM( 'C', 'N', nij, nab, 2*ngm, fact, qgm, 2*ngm, aux, &
            #         2*ngm, 0.0_dp, deeaux, nij )
            #println("sum deeaux = ", sum(deeaux))
            nb = 0
            for ia in 1:Natoms
                if atm2species[ia] != isp
                    continue
                end
                nb = nb + 1
                ijh = 0
                for ih in 1:nh[isp], jh in ih:nh[isp]
                    ijh = ijh + 1
                    Deeq[ih,jh,ia,ispin] = CellVolume * deeaux[ijh,nb]
                    #@printf("%4d %4d %4d %4d %18.10f\n", ih, jh, ia, ispin, Deeq[ih,jh,ia,ispin])
                    if jh > ih
                        Deeq[jh,ih,ia,ispin] = Deeq[ih,jh,ia,ispin]
                    end
                end
            end # Natoms

        end # Nspin
            
    end # Nspecies

    #println("In calc_integ_VQeff: sum Deeq: ", sum(Deeq))

    return

end


function calc_newDeeq!( Ham )

    pspotNL = Ham.pspotNL

    # Early return if no USPP is used
    # XXX probably this should be stored as variable in pspotNL
    ok_paw = any(Ham.pspotNL.are_paw)
    ok_uspp = any(Ham.pspotNL.are_ultrasoft)
    need_augmentation = ok_uspp || ok_paw
    if !need_augmentation
        return
    end

    # Calculate the integral between Qfunc and Veff
    calc_integ_QVeff!( Ham )

    # Some treatment for PAW here ...
    # TODO:

    atm2species = Ham.atoms.atm2species
    Nspin = Ham.electrons.Nspin
    Natoms = Ham.atoms.Natoms
    nh = pspotNL.nh
    Dvan = pspotNL.Dvan
    Deeq = pspotNL.Deeq

    # Add Dvan
    println("sum Deeq after calc_integ_QVeff: ", sum(Deeq))
    println("sum Dvan = ", sum(Dvan))
    #println("Some Deeq")
    for ia in 1:Natoms
        isp = atm2species[ia]
        for ispin in 1:Nspin
            for ih in 1:nh[isp], jh in ih:nh[isp]
                Deeq[ih,jh,ia,ispin] = Deeq[ih,jh,ia,ispin] + Dvan[ih,jh,isp]
                #@printf("%4d %4d %4d %4d %18.10f\n", ih, jh, ia, ispin, Deeq[ih,jh,ia,ispin]*0.5)
                # Factor of 0.5 to match QE result (to 1/Ry?)
                # Unit of D is 1/Ha -> 1/(2Ry)
                Deeq[jh,ih,ia,ispin] = Deeq[ih,jh,ia,ispin]
            end
        end
    end
    
    println("sum Deeq after USPP contrib (before PAW): ", sum(Deeq))

    # Add PAW component
    if ok_paw
        println("Adding ddd_paw to Deeq")
        ddd_paw = Ham.pspotNL.paw.ddd_paw
        for ia in 1:Natoms
            isp = atm2species[ia]
            if !Ham.pspots[isp].is_paw
                continue
            end
            ijh = 0
            for ispin in 1:Nspin, ih in 1:nh[isp], jh in ih:nh[isp]
                ijh = ijh + 1
                Deeq[ih,jh,ia,ispin] += ddd_paw[ijh,ia,ispin]
                Deeq[jh,ih,ia,ispin] = Deeq[ih,jh,ia,ispin] 
            end
        end
    end

    println("sum Deeq before after PAW contrib if any: ", sum(Deeq))

    return

end