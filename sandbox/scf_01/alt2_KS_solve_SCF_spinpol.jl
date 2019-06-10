function alt2_KS_solve_SCF!( Ham::Hamiltonian ;
                        startingwfc=nothing, savewfc=false,
                        betamix = 0.2, NiterMax=100, verbose=true,
                        print_final_ebands=true,
                        print_final_energies=true,
                        check_rhoe=false,
                        use_smearing = false, kT=1e-3,
                        update_psi="LOBPCG", cheby_degree=8,
                        mix_method="simple", MIXDIM=4,
                        print_e_gap=false,
                        etot_conv_thr=1e-6 )

    pw = Ham.pw
    Ngw = pw.gvecw.Ngw
    wk = Ham.pw.gvecw.kpoints.wk
    #
    kpoints = pw.gvecw.kpoints
    Nkpt = kpoints.Nkpt
    #
    Ns = pw.Ns
    Npoints = prod(Ns)
    CellVolume = pw.CellVolume
    dVol = CellVolume/Npoints
    #
    electrons = Ham.electrons
    Nelectrons = electrons.Nelectrons
    Focc = copy(electrons.Focc) # make sure to use the copy
    Nstates = electrons.Nstates
    Nstates_occ = electrons.Nstates_occ
    Nspin = electrons.Nspin
    #
    Nkspin = Nkpt*Nspin

    Nstates_occ = electrons.Nstates_occ

    @assert( Nspin == 2 )

    #
    # Random guess of wave function
    #
    if startingwfc==nothing
        psiks = rand_BlochWavefunc(pw, electrons)
    else
        psiks = startingwfc
    end

    E_GAP_INFO = false
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

    #
    # Calculated electron density from this wave function and update Hamiltonian
    #
    Rhoe = zeros(Float64,Npoints,Nspin)
    
    # introduce different Rhoe_up and Rhoe_dn
    if Nspin == 2
        idx_HOMO = max(round(Int64,Nstates_occ/2),1)
        idx_LUMO = idx_HOMO + 1
        Focc[idx_HOMO,1:Nkpt] = Focc[idx_HOMO,1:Nkpt] .+ 0.5
        Focc[idx_HOMO,Nkpt+1:2*Nkpt] = Focc[idx_HOMO,Nkpt+1:2*Nkpt] .- 0.5
    end

    Rhoe[:,:] = calc_rhoe( Nelectrons, pw, Focc, psiks, Nspin )
    if Nspin == 2
        Rhoe_tot = Rhoe[:,1] + Rhoe[:,2]
        magn_den = Rhoe[:,1] - Rhoe[:,2]
        @printf("\nInitial integ Rhoe up = %18.10f\n", sum(Rhoe[:,1])*dVol)
        @printf("\nInitial integ Rhoe dn = %18.10f\n", sum(Rhoe[:,2])*dVol)
        @printf("\nInitial integ magn_den = %18.10f\n", sum(magn_den)*dVol)
    end

    update!(Ham, Rhoe)

    Etot_old = 0.0

    Rhoe_new = zeros(Float64,Npoints,Nspin)

    diffRhoe = zeros(Nspin)

    evals = zeros(Float64,Nstates,Nkspin)

    ETHR_EVALS_LAST = 1e-6

    ethr = 0.1
    
    if mix_method == "anderson"
        df = zeros(Float64,Npoints*Nspin,MIXDIM)
        dv = zeros(Float64,Npoints*Nspin,MIXDIM)
    
    elseif mix_method in ("rpulay", "rpulay_kerker")
        XX = zeros(Float64,Npoints*Nspin,MIXDIM)
        FF = zeros(Float64,Npoints*Nspin,MIXDIM)
        x_old = zeros(Float64,Npoints*Nspin)
        f_old = zeros(Float64,Npoints*Nspin)

        #=XX = zeros(Float64,Npoints,MIXDIM)
        FF = zeros(Float64,Npoints,MIXDIM)
        x_old = zeros(Float64,Npoints)
        f_old = zeros(Float64,Npoints)

        mgn_XX = zeros(Float64,Npoints,MIXDIM)
        mgn_FF = zeros(Float64,Npoints,MIXDIM)
        mgn_x_old = zeros(Float64,Npoints)
        mgn_f_old = zeros(Float64,Npoints)=#
    end


    @printf("\n")
    @printf("Self-consistent iteration begins ...\n")
    @printf("update_psi = %s\n", update_psi)
    @printf("\n")
    @printf("mix_method = %s\n", mix_method)
    if mix_method in ("rpulay", "rpulay_kerker", "anderson")
        @printf("MIXDIM = %d\n", MIXDIM)
    end
    @printf("Density mixing with betamix = %10.5f\n", betamix)
    if use_smearing
        @printf("Smearing = %f\n", kT)
    end
    println("") # blank line before SCF iteration info

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( Ham.atoms )

    # calculate PspCore energy
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )

    CONVERGED = 0

    E_fermiSpin = zeros(Nspin)

    for iter = 1:NiterMax

        # determine convergence criteria for diagonalization
        if iter == 1
            ethr = 0.1
        elseif iter == 2
            ethr = 0.01
        else
            ethr = ethr/5.0
            ethr = max( ethr, ETHR_EVALS_LAST )
        end

        if update_psi == "LOBPCG"

            evals =
            diag_LOBPCG!( Ham, psiks, verbose=false, verbose_last=false,
                          Nstates_conv=Nstates_occ )

        elseif update_psi == "davidson"

            evals =
            diag_davidson!( Ham, psiks, verbose=false, verbose_last=false,
                            Nstates_conv=Nstates_occ )                

        elseif update_psi == "PCG"

            evals =
            diag_Emin_PCG!( Ham, psiks, verbose=false, verbose_last=false,
                            Nstates_conv=Nstates_occ )

        elseif update_psi == "CheFSI"
            # evals will be calculated later
            diag_CheFSI!( Ham, psiks, cheby_degree )

        else
            error( @sprintf("Unknown method for update_psi = %s\n", update_psi) )
        end

        if E_GAP_INFO && verbose && print_e_gap
            @printf("E gap = %18.10f\n", minimum(evals[idx_LUMO,:] - evals[idx_HOMO,:]))
        end

        if use_smearing
            Focc, E_fermi = calc_Focc( Nelectrons, wk, kT, evals, Nspin )
            Entropy = calc_entropy( wk, kT, evals, E_fermi, Nspin )
            Ham.electrons.Focc = copy(Focc)
        end

        Rhoe_new[:,:] = calc_rhoe( Nelectrons, pw, Focc, psiks, Nspin )
        for ispin = 1:Nspin
            diffRhoe[ispin] = norm(Rhoe_new[:,ispin] - Rhoe[:,ispin])
        end

        # check norm of
        if check_rhoe
            integRhoe = sum(Rhoe_new)*dVol
            @printf("Before: integRhoe_new = %18.10f\n", integRhoe)
            Rhoe_new = Nelectrons/integRhoe * Rhoe_new
            integRhoe = sum(Rhoe_new)*dVol
            @printf("After renormalize Rhoe_new: = %18.10f\n", integRhoe)
        end

        if mix_method == "simple"
            #for ispin = 1:Nspin
            #    Rhoe[:,ispin] = betamix*Rhoe_new[:,ispin] + (1-betamix)*Rhoe[:,ispin]
            #end

            # FIXME: only works for Nspin == 2
            Rhoe_tot_new = Rhoe_new[:,1] + Rhoe_new[:,2]
            magn_den_new = Rhoe_new[:,1] - Rhoe_new[:,2]
    
            Rhoe_tot = betamix*Rhoe_tot_new + (1-betamix)*Rhoe_tot
            magn_den = betamix*magn_den_new + (1-betamix)*magn_den
    
            Rhoe[:,1] = 0.5*(Rhoe_tot + magn_den)
            Rhoe[:,2] = 0.5*(Rhoe_tot - magn_den)
    

        elseif mix_method == "simple_kerker"
            for ispin = 1:Nspin
                Rhoe[:,ispin] = Rhoe[:,ispin] + betamix*precKerker(pw, Rhoe_new[:,ispin] - Rhoe[:,ispin])
            end

        elseif mix_method == "rpulay"
            
            Rhoe = reshape( mix_rpulay!(
                reshape(Rhoe,(Npoints*Nspin)),
                reshape(Rhoe_new,(Npoints*Nspin)), betamix, XX, FF, iter, MIXDIM, x_old, f_old
                ), (Npoints,Nspin) )

            if Nspin == 2
                magn_den = Rhoe[:,1] - Rhoe[:,2]
            end
            

            #=
            Rhoe_tot_new = Rhoe_new[:,1] + Rhoe_new[:,2]
            magn_den_new = Rhoe_new[:,1] - Rhoe_new[:,2]
        
            Rhoe_tot = mix_rpulay!(
                Rhoe_tot, Rhoe_tot_new, betamix,
                XX, FF, iter, MIXDIM, x_old, f_old )
            
            magn_den = mix_rpulay!(
                magn_den, magn_den_new, betamix,
                mgn_XX, mgn_FF, iter, MIXDIM, mgn_x_old, mgn_f_old )
        
            Rhoe[:,1] = 0.5*(Rhoe_tot + magn_den)
            Rhoe[:,2] = 0.5*(Rhoe_tot - magn_den)
            =#
            
        elseif mix_method == "rpulay_kerker"
        
            Rhoe = reshape( mix_rpulay_kerker!( pw,
                reshape(Rhoe,(Npoints*Nspin)),
                reshape(Rhoe_new,(Npoints*Nspin)), betamix, XX, FF, iter, MIXDIM, x_old, f_old
                ), (Npoints,Nspin) )
            
            if Nspin == 2
                magn_den = Rhoe[:,1] - Rhoe[:,2]
            end
        
        elseif mix_method == "anderson"
            Rhoe[:,:] = mix_anderson!( Nspin, Rhoe, Rhoe_new, betamix, df, dv, iter, MIXDIM )
        
        else
            error(@sprintf("Unknown mix_method = %s\n", mix_method))

        end

        for rho in Rhoe
            if rho < eps()
                rho = 0.0
            end
        end

        # renormalize
        if check_rhoe
            integRhoe = sum(Rhoe)*dVol
            @printf("After mixing: integRhoe = %18.10f\n", integRhoe)
            Rhoe = Nelectrons/integRhoe * Rhoe
            integRhoe = sum(Rhoe)*dVol
            @printf("After renormalize Rhoe: = %18.10f\n", integRhoe)
        end

        update!( Ham, Rhoe )

        # Calculate energies
        Ham.energies = calc_energies( Ham, psiks )
        if use_smearing
            Ham.energies.mTS = Entropy
        end
        Etot = sum(Ham.energies)
        diffE = abs( Etot - Etot_old )

        if verbose
            if Nspin == 1
                @printf("SCF: %8d %18.10f %18.10e %18.10e\n",
                        iter, Etot, diffE, diffRhoe[1] )
            else
                @printf("SCF: %8d %18.10f %18.10e %18.10e %18.10e\n",
                        iter, Etot, diffE, diffRhoe[1], diffRhoe[2] )
                magn_den = Rhoe[:,1] - Rhoe[:,2]
                @printf("integ Rhoe spin up = %18.10f\n", sum(Rhoe[:,1])*dVol) 
                @printf("integ Rhoe spin dn = %18.10f\n", sum(Rhoe[:,2])*dVol) 
                @printf("integ magn_den = %18.10f\n", sum(magn_den)*dVol) 
            end
        
        end

        if diffE < etot_conv_thr
            CONVERGED = CONVERGED + 1
        else  # reset CONVERGED
            CONVERGED = 0
        end

        if CONVERGED >= 2
            if verbose
                @printf("SCF is converged: iter: %d , diffE = %10.7e\n", iter, diffE)
            end
            break
        end
        #
        Etot_old = Etot

        flush(stdout)
    end

    # Eigenvalues are not calculated if using CheFSI.
    # We calculate them here.
    if update_psi == "CheFSI"
        for ispin = 1:Nspin
        for ik = 1:Nkpt
            Ham.ik = ik
            Ham.ispin = ispin
            ikspin = ik + (ispin - 1)*Nkpt
            Hr = psiks[ikspin]' * op_H( Ham, psiks[ikspin] )
            evals[:,ikspin] = eigvals(Hermitian(Hr))
        end
        end
    end

    Ham.electrons.ebands = evals

    if verbose && print_final_ebands
        @printf("\n")
        @printf("----------------------------\n")
        @printf("Final Kohn-Sham eigenvalues:\n")
        @printf("----------------------------\n")
        @printf("\n")
        print_ebands(Ham.electrons, Ham.pw.gvecw.kpoints)
    end

    if verbose && print_final_energies
        @printf("\n")
        @printf("-------------------------\n")
        @printf("Final Kohn-Sham energies:\n")
        @printf("-------------------------\n")
        @printf("\n")
        println(Ham.energies, use_smearing=use_smearing)
    end

    if savewfc
        for ikspin = 1:Nkpt*Nspin
            wfc_file = open("WFC_ikspin_"*string(ikspin)*".data","w")
            write( wfc_file, psiks[ikspin] )
            close( wfc_file )
        end
    end

    return

end
