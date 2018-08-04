function KS_solve_SCF_smearing!( Ham::Hamiltonian ;
                       startingwfc=nothing, savewfc=false,
                       betamix = 0.5, NiterMax=100, verbose=false,
                       check_rhoe_after_mix=false,
                       update_psi="LOBPCG", cheby_degree=8,
                       mix_method="simple",
                       ETOT_CONV_THR=1e-6 )

    pw = Ham.pw
    Ngw = pw.gvecw.Ngw
    wk = Ham.pw.gvecw.kpoints.wk
    #
    kpoints = pw.gvecw.kpoints
    Nkpt = kpoints.Nkpt
    #
    Ns = pw.Ns
    Npoints = prod(Ns)
    dVol = pw.CellVolume/Npoints
    #
    electrons = Ham.electrons
    Nelectrons = electrons.Nelectrons
    Focc_orig = electrons.Focc
    Nstates = electrons.Nstates
    Nspin = electrons.Nspin
    #
    Nkspin = Nkpt*Nspin

    Nstates_occ = electrons.Nstates_occ

    psiks = Array{Array{ComplexF64,2},1}(undef,Nkspin)

    #
    # Random guess of wave function
    #
    if startingwfc==nothing
        Base.seed!(1234)
        for ispin = 1:Nspin
        for ik = 1:Nkpt
            ikspin = ik + (ispin - 1)*Nkpt
            psi = rand(Ngw[ik],Nstates) + im*rand(Ngw[ik],Nstates)
            psiks[ikspin] = ortho_gram_schmidt(psi)
        end
        end
    else
        psiks = startingwfc
    end

    #
    # Calculated electron density from this wave function and update Hamiltonian
    #
    Rhoe = zeros(Float64,Npoints,Nspin)
    for ispin = 1:Nspin
        idxset = (Nkpt*(ispin-1)+1):(Nkpt*ispin)
        Rhoe[:,ispin] = calc_rhoe( pw, Focc_orig[:,idxset], psiks[idxset] )
    end
    update!(Ham, Rhoe)

    Etot_old = 0.0

    Rhoe_new = zeros(Float64,Npoints,Nspin)

    diffRhoe = zeros(Nspin)

    evals = zeros(Float64,Nstates,Nkspin)

    ETHR_EVALS_LAST = 1e-6

    ethr = 0.1

    # For Anderson mixing
    MIXDIM = 4
    if mix_method == "anderson"
        df = zeros(Float64,Npoints*Nspin,MIXDIM)
        dv = zeros(Float64,Npoints*Nspin,MIXDIM)
    end

    E_GAP_INFO = false
    Nstates_occ = electrons.Nstates_occ
    if Nstates_occ < Nstates
        E_GAP_INFO = true
        if Nspin == 2
            idx_HOMO = max(round(Int64,Nstates_occ/2),1)
            idx_LUMO = idx_HOMO + 1
        else
            idx_HOMO = Nstates_occ
            idx_LUMO = idx_HOMO + 1
        end
    end
    
    Focc = zeros(Nstates,Nkspin)
    kT = 0.01

    @printf("\n")
    @printf("Self-consistent iteration begins ...\n")
    @printf("\n")
    if mix_method == "anderson"
        @printf("Using Anderson mixing\n")
    else
        @printf("Using simple mixing\n")
    end
    @printf("Density mixing with betamix = %10.5f\n", betamix)
    @printf("Smearing = %f\n", kT)
    @printf("\n")

    CONVERGED = 0

    for iter = 1:NiterMax

        println("")

        if update_psi == "LOBPCG"
            for ispin = 1:Nspin
            for ik = 1:Nkpt
                Ham.ik = ik
                Ham.ispin = ispin
                ikspin = ik + (ispin - 1)*Nkpt
                #
                evals[:,ikspin], psiks[ikspin] =
                diag_LOBPCG( Ham, psiks[ikspin], verbose=false, verbose_last=false,
                             Nstates_conv = Nstates_occ )
                #
            end
            end

        elseif update_psi == "PCG"
            
            # determined convergence criteria for diagonalization
            if iter == 1
                ethr = 0.1
            elseif iter == 2
                ethr = 0.01
            else
                ethr = ethr/5.0
                ethr = max( ethr, ETHR_EVALS_LAST )
            end

            for ispin = 1:Nspin
            for ik = 1:Nkpt
                Ham.ik = ik
                Ham.ispin = ispin
                ikspin = ik + (ispin - 1)*Nkpt
                evals[:,ikspin], psiks[ikspin] =
                diag_Emin_PCG( Ham, psiks[ikspin], TOL_EBANDS=ethr, verbose=false )
            end
        end

        elseif update_psi == "CheFSI"
            
            for ik = 1:Nkpt
                Ham.ik = ik
                Ham.ispin = ispin
                ikspin = ik + (ispin - 1)*Nkpt
                ub, lb = get_ub_lb_lanczos( Ham, Nstates*2 )
                psiks[ikspin] = chebyfilt( Ham, psiks[ikspin], cheby_degree, lb, ub)
                psiks[ikspin] = ortho_gram_schmidt( psiks[ik] )
            end

        else

            @printf("Unknown method for update_psi = %s\n", update_psi)
            exit()

        end

        if E_GAP_INFO
            println("E gap = ", minimum(evals[idx_LUMO,:] - evals[idx_HOMO,:]))
        end

        Focc, E_fermi = calc_Focc( evals, wk, Nelectrons, kT, Nspin=Nspin )
        Entropy = calc_entropy( Focc, wk, kT, Nspin=Nspin )

        for ispin = 1:Nspin
            idxset = (Nkpt*(ispin-1)+1):(Nkpt*ispin)
            Rhoe_new[:,ispin] = calc_rhoe( pw, Focc[:,idxset], psiks[idxset] )
            diffRhoe[ispin] = norm(Rhoe_new[:,ispin] - Rhoe[:,ispin])
        end

        if mix_method == "simple"
            for ispin = 1:Nspin
                Rhoe[:,ispin] = betamix*Rhoe_new[:,ispin] + (1-betamix)*Rhoe[:,ispin]
            end
        elseif mix_method == "anderson"
            # FIXME: df and dv is not modified when we call it by df[:,:] or dv[:,:]
            #Rhoe[:,:] = andersonmix!( Rhoe, Rhoe_new, betamix, df, dv, iter, MIXDIM )
            Rhoe[:,:] = mix_anderson!( Nspin, Rhoe, Rhoe_new, betamix, df, dv, iter, MIXDIM )
        else
            @printf("ERROR: Unknown mix_method = %s\n", mix_method)
            exit()
        end

        for ispin = 1:Nspin
            for ip = 1:Npoints
                if Rhoe[ip,ispin] < 1e-12
                    Rhoe[ip,ispin] = 1e-12
                end
            end
        end

        # renormalize
        if check_rhoe_after_mix
            integRhoe = sum(Rhoe)*dVol
            @printf("After mixing: integRhoe = %18.10f\n", integRhoe)
            Rhoe = Nelectrons/integRhoe * Rhoe
            integRhoe = sum(Rhoe)*dVol
            @printf("After renormalize Rhoe: = %18.10f\n", integRhoe)
        end

        update!( Ham, Rhoe )

        # Calculate energies
        Ham.energies = calc_energies( Ham, psiks )
        Etot = sum(Ham.energies) + Entropy
        diffE = abs( Etot - Etot_old )

        if Nspin == 1
            @printf("SCF: %8d %18.10f %18.10e %18.10e\n",
                    iter, Etot, diffE, diffRhoe[1] )
        else
            @printf("SCF: %8d %18.10f %18.10e %18.10e %18.10e\n",
                    iter, Etot, diffE, diffRhoe[1], diffRhoe[2] )
        end
        @printf("Entropy (-TS) = %18.10f\n", Entropy)

        if diffE < ETOT_CONV_THR
            CONVERGED = CONVERGED + 1
        else  # reset CONVERGED
            CONVERGED = 0
        end

        if CONVERGED >= 2
            @printf("SCF is converged: iter: %d , diffE = %10.7e\n", iter, diffE)
            break
        end
        #
        Etot_old = Etot
    end

    # Eigenvalues are not calculated if using CheFSI.
    # We calculate them here.
    if update_psi == "CheFSI"
        for ispin = 1:Nspin
        for ik = 1:Nkpt
            Ham.ik = ik
            Ham.ispin = ispin
            ikspin = ik + (ispin - 1)*Nkpt
            Hr = psik[ikspin]' * op_H( Ham, psiks[ikspin] )
            evals[:,ikspin] = real(eigvals(Hr))
        end
        end
    end

    Ham.electrons.ebands = evals

    if savewfc
        for ikspin = 1:Nkpt*Nspin
            wfc_file = open("WFC_ikspin_"*string(ikspin)*".data","w")
            write( wfc_file, psiks[ikspin] )
            close( wfc_file )
        end
    end

    return

end
