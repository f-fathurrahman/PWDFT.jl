"""
Solves Kohn-Sham problem using direct energy minimization as described
by Ismail-Beigi and Arias.
"""
function KS_solve_Emin_PCG!( Ham::Hamiltonian;
                             startingwfc=:random, savewfc=false,
                             startingrhoe=:gaussian,
                             skip_initial_diag=false,
                             α_t=3e-5, NiterMax=200, verbose=true,
                             print_final_ebands=false, print_final_energies=true,
                             i_cg_beta=2, etot_conv_thr=1e-6 )

    pw = Ham.pw
    electrons = Ham.electrons
    
    Focc = electrons.Focc
    Nstates = electrons.Nstates
    Nelectrons = electrons.Nelectrons
    
    Ns = pw.Ns
    Npoints = prod(Ns)
    CellVolume = pw.CellVolume
    Ngw = pw.gvecw.Ngw
    Ngwx = pw.gvecw.Ngwx
    Nkpt = pw.gvecw.kpoints.Nkpt
    
    Nspin = electrons.Nspin
    Nkspin = Nkpt*Nspin

    #
    # Initial wave function
    #
    if startingwfc == :read
        psiks = read_psiks( Ham )
    else
        # generate random BlochWavefunc
        psiks = rand_BlochWavefunc( Ham )
    end

    if Ham.sym_info.Nsyms > 1
        rhoe_symmetrizer = RhoeSymmetrizer( Ham )
    end

    #
    # Calculated electron density from this wave function and
    # update Hamiltonian (calculate Hartree and XC potential).
    #
    Rhoe = zeros(Float64,Npoints,Nspin)

    if startingrhoe == :gaussian
        @assert Nspin == 1
        Rhoe[:,1] = guess_rhoe( Ham )
        #Rhoe = guess_rhoe_atomic( Ham ) # when smearing is ready this should be preferred
    else
        calc_rhoe!( Ham, psiks, Rhoe )
    end
    
    # Symmetrize Rhoe is needed
    if Ham.sym_info.Nsyms > 1
        symmetrize_rhoe!( Ham, rhoe_symmetrizer, Rhoe )
    end

    update!(Ham, Rhoe)

    evals = zeros(Nstates,Nkspin)
    if !skip_initial_diag
        evals =
        diag_LOBPCG!( Ham, psiks, verbose=false, verbose_last=false, NiterMax=10 )
    end

    #
    # Variables for PCG
    #
    g      = zeros_BlochWavefunc( Ham )
    d      = zeros_BlochWavefunc( Ham )
    g_old  = zeros_BlochWavefunc( Ham )
    d_old  = zeros_BlochWavefunc( Ham )
    Kg     = zeros_BlochWavefunc( Ham )
    Kg_old = zeros_BlochWavefunc( Ham )
    psic   = zeros_BlochWavefunc( Ham )
    gt     = zeros_BlochWavefunc( Ham )
    
    β = zeros(Nkspin)
    α = zeros(Nkspin)

    Etot_old = 0.0

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( Ham.atoms )

    # calculate PspCore energy
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )

    # Calculate energy at this psi
    energies = calc_energies(Ham, psiks)
    Ham.energies = energies
    Etot = sum(energies)

    CONVERGED = 0

    if verbose
        @printf("\n")
        @printf("Minimizing Kohn-Sham energy using PCG\n")
        @printf("-------------------------------------\n")
        @printf("NiterMax  = %d\n", NiterMax)
        @printf("α_t       = %e\n", α_t)
        @printf("conv_thr  = %e\n", etot_conv_thr)
        if i_cg_beta == 1
            @printf("Using Fletcher-Reeves formula for CG_BETA\n")
        elseif i_cg_beta == 2
            @printf("Using Polak-Ribiere formula for CG_BETA\n")
        elseif i_cg_beta == 3
            @printf("Using Hestenes-Stiefeld formula for CG_BETA\n")
        else
            @printf("Using Dai-Yuan formula for CG_BETA\n")
        end
        @printf("\n")
    end


    for iter = 1:NiterMax

        for ispin = 1:Nspin
        for ik = 1:Nkpt

            Ham.ik = ik
            Ham.ispin = ispin
            ikspin = ik + (ispin - 1)*Nkpt

            g[ikspin] = calc_grad( Ham, psiks[ikspin] )
            Kg[ikspin] = Kprec( Ham.ik, pw, g[ikspin] )

            # XXX: define function trace for real(sum(conj(...)))
            if iter != 1
                if i_cg_beta == 1
                    β[ikspin] =
                    real(sum(conj(g[ikspin]).*Kg[ikspin]))/real(sum(conj(g_old[ikspin]).*Kg_old[ikspin]))
                elseif i_cg_beta == 2
                    β[ikspin] =
                    real(sum(conj(g[ikspin]-g_old[ikspin]).*Kg[ikspin]))/real(sum(conj(g_old[ikspin]).*Kg_old[ikspin]))
                elseif i_cg_beta == 3
                    β[ikspin] =
                    real(sum(conj(g[ikspin]-g_old[ikspin]).*Kg[ikspin]))/real(sum(conj(g[ikspin]-g_old[ikspin]).*d[ikspin]))
                else
                    β[ikspin] =
                    real(sum(conj(g[ikspin]).*Kg[ikspin]))/real(sum((g[ikspin]-g_old[ikspin]).*conj(d_old[ikspin])))
                end
            end
            if β[ikspin] < 0.0
                β[ikspin] = 0.0
            end

            d[ikspin] = -Kg[ikspin] + β[ikspin] * d_old[ikspin]

            psic[ikspin] = ortho_sqrt(psiks[ikspin] + α_t*d[ikspin])
        end # ik
        end # ispin
        
        calc_rhoe!( Ham, psiks, Rhoe )
        # Symmetrize Rhoe is needed
        if Ham.sym_info.Nsyms > 1
            symmetrize_rhoe!( Ham, rhoe_symmetrizer, Rhoe )
        end

        update!(Ham, Rhoe)

        for ispin = 1:Nspin
        for ik = 1:Nkpt
            Ham.ik = ik
            Ham.ispin = ispin
            ikspin = ik + (ispin - 1)*Nkpt
            gt[ikspin] = calc_grad(Ham, psic[ikspin])

            denum = real(sum(conj(g[ikspin]-gt[ikspin]).*d[ikspin]))
            if denum != 0.0
                α[ikspin] = abs( α_t*real(sum(conj(g[ikspin]).*d[ikspin]))/denum )
            else
                α[ikspin] = 0.0
            end

            # Update wavefunction
            psiks[ikspin] = psiks[ikspin] + α[ikspin]*d[ikspin]

            # Update potentials
            psiks[ikspin] = ortho_sqrt(psiks[ikspin])
        end
        end

        calc_rhoe!( Ham, psiks, Rhoe )
        # Symmetrize Rhoe is needed
        if Ham.sym_info.Nsyms > 1
            symmetrize_rhoe!( Ham, rhoe_symmetrizer, Rhoe )
        end
        update!(Ham, Rhoe)

        Ham.energies = calc_energies( Ham, psiks )
        Etot = sum(Ham.energies)
        diffE = abs(Etot-Etot_old)

        if verbose
            @printf("CG step %8d = %18.10f %10.7e\n", iter, Etot, diffE)
        end
        
        if diffE < etot_conv_thr
            CONVERGED = CONVERGED + 1
        else
            CONVERGED = 0
        end

        if CONVERGED >= 2
            if verbose
                @printf("CONVERGENCE ACHIEVED\n")
            end
            break
        end

        g_old = copy(g)
        d_old = copy(d)
        Kg_old = copy(Kg)
        Etot_old = Etot

        flush(stdout)
    end

    # Calculate eigenvalues
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        psiks[ikspin] = ortho_sqrt(psiks[ikspin])
        Hr = Hermitian(psiks[ikspin]' * op_H(Ham, psiks[ikspin]))
        evals, evecs = eigen(Hr)
        Ham.electrons.ebands[:,ik] = evals
        psiks[ikspin] = psiks[ikspin]*evecs
    end
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
        println(Ham.energies)
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
