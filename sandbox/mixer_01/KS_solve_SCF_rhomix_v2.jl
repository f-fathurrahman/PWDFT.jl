function KS_solve_SCF_rhomix_v2!(
    Ham::Hamiltonian ;
    startingwfc=:random,
    savewfc=false,
    startingrhoe=:gaussian,
    betamix=0.2,
    NiterMax=100,
    verbose=true,
    print_final_ebands=false,
    print_final_energies=true,
    print_integ_rhoe=false,
    check_rhoe=false,
    use_smearing=false,
    kT=1e-3,
    update_psi="LOBPCG",
    cheby_degree=8,
    mix_method="simple",
    mixdim=5,
    print_e_gap=false,
    etot_conv_thr=1e-6,
    ethr_evals_last=1e-5,
    starting_magnetization=nothing
)

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

    #
    # Initial wave function
    #
    if startingwfc == :read
        psiks = read_psiks( Ham )
    else
        # generate random BlochWavefunc
        psiks = rand_BlochWavefunc( Ham )
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

    # Workspace for Rhoe symmetrization is initialized here
    if Ham.sym_info.Nsyms > 1
        rhoe_symmetrizer = RhoeSymmetrizer( Ham )
    end

    #
    # Calculated electron density from this wave function and update Hamiltonian
    #
    Rhoe = zeros(Float64,Npoints,Nspin)

    if startingrhoe == :gaussian && startingwfc == :random
        if Nspin == 1
            Rhoe[:,1] = guess_rhoe( Ham )
        else
            Rhoe = guess_rhoe_atomic( Ham, starting_magnetization=starting_magnetization )
        end
    else
        Rhoe = calc_rhoe( Nelectrons, pw, Focc, psiks, Nspin )
    end
    # Symmetrize Rhoe if needed
    if Ham.sym_info.Nsyms > 1
        symmetrize_rhoe!( Ham, rhoe_symmetrizer, Rhoe )
    end

    if Nspin == 2
        magn_den = zeros(Npoints)
        for ip in 1:Npoints
            magn_den[ip] = Rhoe[ip,1] - Rhoe[ip,2]
        end
    end

    if Nspin == 2 && verbose
        @printf("Initial integ Rhoe up  = %18.10f\n", sum(Rhoe[:,1])*dVol)
        @printf("Initial integ Rhoe dn  = %18.10f\n", sum(Rhoe[:,2])*dVol)
        @printf("Initial integ magn_den = %18.10f\n", sum(magn_den)*dVol)
    end

    update!(Ham, Rhoe)

    Etot_old = 0.0

    Rhoe_new = zeros(Float64,Npoints,Nspin)

    diffRhoe = zeros(Nspin)

    evals = zeros(Float64,Nstates,Nkspin)

    ethr = 0.1
    
    if mix_method in ("anderson", "broyden")
        df = zeros(Float64,Npoints*Nspin, mixdim)
        dv = zeros(Float64,Npoints*Nspin, mixdim)

    elseif mix_method == "linear_adaptive"
        betav = betamix*ones(Float64, Npoints*Nspin)
        df = zeros(Float64, Npoints*Nspin)
    
    elseif mix_method in ("rpulay", "ppulay", "pulay")
        XX = zeros(Float64,Npoints*Nspin, mixdim)
        FF = zeros(Float64,Npoints*Nspin, mixdim)
        x_old = zeros(Float64,Npoints,Nspin)
        f_old = zeros(Float64,Npoints,Nspin)
    end


    @printf("\n")
    @printf("Self-consistent iteration begins ...\n")
    @printf("update_psi = %s\n", update_psi)
    @printf("mix_method = %s\n", mix_method)
    if mix_method in ("rpulay", "anderson", "ppulay", "broyden")
        @printf("mixdim = %d\n", mixdim)
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

    Nconverges = 0

    E_fermi = 0.0

    if verbose
        if Nspin == 1
            @printf("-------------------------------------------------------\n")
            @printf("       iter            E            ΔE           Δρ\n")
            @printf("-------------------------------------------------------\n")
        else
            @printf("----------------------------------------------------------------------\n")
            @printf("       iter            E            ΔE                  Δρ\n")
            @printf("----------------------------------------------------------------------\n")
        end
    end

    for iter = 1:NiterMax

        # determine convergence criteria for diagonalization
        if iter == 1
            ethr = 0.1
        elseif iter == 2
            ethr = 0.01
        else
            ethr = ethr/5.0
            ethr = max( ethr, ethr_evals_last )
        end

        if update_psi == "LOBPCG"

            evals =
            diag_LOBPCG!( Ham, psiks, verbose=false, verbose_last=false, tol=ethr,
                          Nstates_conv=Nstates_occ )

        elseif update_psi == "davidson"

            evals =
            diag_davidson!( Ham, psiks, verbose=false, verbose_last=false, tol=ethr,
                            Nstates_conv=Nstates_occ )                

        elseif update_psi == "PCG"

            evals =
            diag_Emin_PCG!( Ham, psiks, verbose=false, verbose_last=false, tol=ethr,
                            Nstates_conv=Nstates_occ )

        elseif update_psi == "CheFSI"

            evals =
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
        if Ham.sym_info.Nsyms > 1
            symmetrize_rhoe!( Ham, rhoe_symmetrizer, Rhoe_new )
        end

        for ispin = 1:Nspin
            diffRhoe[ispin] = sum(abs.(Rhoe_new[:,ispin] - Rhoe[:,ispin]))/Npoints
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

            Rhoe = betamix*Rhoe_new + (1-betamix)*Rhoe

        elseif mix_method == "linear_adaptive"

            mix_adaptive!( Rhoe, Rhoe_new, betamix, betav, df )

        elseif mix_method == "broyden"

            mix_broyden!( Rhoe, Rhoe_new, betamix, iter, mixdim, df, dv )

        elseif mix_method == "pulay"
        
            mix_pulay!( Rhoe, Rhoe_new, betamix, XX, FF, iter, mixdim, x_old, f_old )

        elseif mix_method == "rpulay"
            
            mix_rpulay!( Rhoe, Rhoe_new, betamix, XX, FF, iter, mixdim, x_old, f_old )
            # result is in Rhoe

        elseif mix_method == "ppulay"
            
            #XXX We fix the period to be 3 here

            mix_ppulay!( Rhoe, Rhoe_new, betamix, XX, FF, iter, mixdim, 3, x_old, f_old )

        
        elseif mix_method == "anderson"

            mix_anderson!( Rhoe, Rhoe_new, betamix, df, dv, iter, mixdim )
        
        else
            error(@sprintf("Unknown mix_method = %s\n", mix_method))

        end

        for rho in Rhoe
            if rho < eps()
                rho = 0.0
            end
        end

        if Nspin == 2
            magn_den = Rhoe[:,1] - Rhoe[:,2]
        end

        # renormalize
        if check_rhoe
            integRhoe = sum(Rhoe)*dVol
            #@printf("After mixing: integRhoe = %18.10f\n", integRhoe)
            Rhoe = Nelectrons/integRhoe * Rhoe
            integRhoe = sum(Rhoe)*dVol
            #@printf("After renormalize Rhoe: = %18.10f\n", integRhoe)
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
                @printf("SCF: %5d %18.10f %12.5e %12.5e\n",
                         iter, Etot, diffE, diffRhoe[1] )
                if print_integ_rhoe
                    @printf("integ Rhoe = %18.10f\n", sum(Rhoe)*dVol)
                end
            else
                @printf("SCF: %5d %18.10f %12.5e [%12.5e,%12.5e]\n",
                         iter, Etot, diffE, diffRhoe[1], diffRhoe[2] )
                if print_integ_rhoe
                    magn_den = Rhoe[:,1] - Rhoe[:,2]
                    @printf("integ Rhoe spin up = %18.10f\n", sum(Rhoe[:,1])*dVol) 
                    @printf("integ Rhoe spin dn = %18.10f\n", sum(Rhoe[:,2])*dVol) 
                    @printf("integ magn_den = %18.10f\n", sum(magn_den)*dVol) 
                end
            end     
        end

        if diffE < etot_conv_thr
            Nconverges = Nconverges + 1
        else
            Nconverges = 0
        end

        if Nconverges >= 2
            if verbose
                @printf("\nSCF is converged in iter: %d\n", iter)
            end
            break
        end
        #
        Etot_old = Etot

        flush(stdout)
    end

    Ham.electrons.ebands = evals

    if use_smearing && verbose
        @printf("\nFermi energy = %18.10f Ha = %18.10f eV\n", E_fermi, E_fermi*2*Ry2eV)
    end

    if Nspin == 2 && verbose
        @printf("\n")
        @printf("Final integ Rhoe up  = %18.10f\n", sum(Rhoe[:,1])*dVol)
        @printf("Final integ Rhoe dn  = %18.10f\n", sum(Rhoe[:,2])*dVol)
        @printf("Final integ magn_den = %18.10f\n", sum(magn_den)*dVol)
    end

    if verbose && print_final_ebands
        @printf("\n")
        @printf("----------------------------\n")
        @printf("Final Kohn-Sham eigenvalues:\n")
        @printf("----------------------------\n")
        @printf("\n")
        print_ebands(Ham.electrons, Ham.pw.gvecw.kpoints, unit="eV")
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
