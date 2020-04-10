import PWDFT: KS_solve_Emin_PCG!

function KS_solve_Emin_PCG!( Ham::CuHamiltonian;
                             startingwfc=:random,
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
        error("startingwfc == :read is not yet implemented")
        #psiks = read_psiks( Ham )
    else
        # generate random BlochWavefunc
        psiks = rand_CuBlochWavefunc( Ham )
    end

    #
    # Calculated electron density from this wave function and
    # update Hamiltonian (calculate Hartree and XC potential).
    #
    Rhoe = CuArrays.zeros(Float64,Npoints,Nspin)

    if startingrhoe == :gaussian
        error("startingrhoe == :gaussian is not yet implemented")
        @assert Nspin == 1
        Rhoe[:,1] = guess_rhoe( Ham )
        #Rhoe = guess_rhoe_atomic( Ham ) # when smearing is ready this should be preferred
    else
        calc_rhoe!( Ham, psiks, Rhoe )
    end

    update!(Ham, Rhoe)

    evals = CuArrays.zeros(Nstates,Nkspin)
    if !skip_initial_diag
        error("skip_initial_diag == false is not yet implemented")
        evals =
        diag_LOBPCG!( Ham, psiks, verbose=false, verbose_last=false, NiterMax=10 )
    end

    # warm-up CuArrays.CUSOLVER
    vv = CuArrays.rand(ComplexF64,1,1)
    _, _ = eigen_heevd_gpu( vv )

    #
    # Variables for PCG
    #
    g      = zeros_CuBlochWavefunc( Ham )
    d      = copy(g)
    g_old  = copy(g)
    d_old  = copy(g)
    Kg     = copy(g)
    Kg_old = copy(g)
    psic   = copy(g)
    gt     = copy(g)
    
    β = zeros(Nkspin)
    α = zeros(Nkspin)

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( Ham.atoms )

    # Calculate energy at this psi
    energies = calc_energies(Ham, psiks)
    Ham.energies = energies
    Etot_old = sum(energies)

    Nconverges = 0

    if verbose
        @printf("\n")
        @printf("Minimizing Kohn-Sham energy using PCG (CUDA version)\n")
        @printf("----------------------------------------------------\n")
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

    CuArrays.memory_status()
    println()

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
                    if abs(β[ikspin] - 1.0) < 0.2
                        β[ikspin] = 0.0
                    end
                elseif i_cg_beta == 2
                    β[ikspin] =
                    real(sum(conj(g[ikspin]-g_old[ikspin]).*Kg[ikspin]))/real(sum(conj(g_old[ikspin]).*Kg_old[ikspin]))
                elseif i_cg_beta == 3
                    β[ikspin] =
                    real(sum(conj(g[ikspin]-g_old[ikspin]).*Kg[ikspin]))/real(sum(conj(g[ikspin]-g_old[ikspin]).*d[ikspin]))
                else
                    β[ikspin] =
                    real(sum(conj(g[ikspin]).*Kg[ikspin]))/real(sum((g[ikspin]-g_old[ikspin]).*conj(d_old[ikspin])))
                    if abs(β[ikspin] - 1.0) < 0.2
                        β[ikspin] = 0.0
                    end                    
                end
            end
            if β[ikspin] < 0.0
                β[ikspin] = 0.0
            end

            d[ikspin] = -Kg[ikspin] + β[ikspin] * d_old[ikspin]

            psic[ikspin] = ortho_gram_schmidt( psiks[ikspin] + α_t*d[ikspin] )
        end # ik
        end # ispin
        
        calc_rhoe!( Ham, psic, Rhoe )

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
            #psiks[ikspin] = ortho_gram_schmidt(psiks[ikspin])
            ortho_gram_schmidt!( psiks[ikspin] )
        end
        end

        # Update Rhoe and potentials
        calc_rhoe!( Ham, psiks, Rhoe )

        update!(Ham, Rhoe)

        Ham.energies = calc_energies( Ham, psiks )
        Etot = sum(Ham.energies)
        diffE = abs(Etot-Etot_old)

        if verbose
            @printf("Emin_PCG step %8d = %18.10f %10.7e\n", iter, Etot, diffE)
        end
        
        if diffE < etot_conv_thr
            Nconverges = Nconverges + 1
        else
            Nconverges = 0
        end

        if Nconverges >= 2
            if verbose
                @printf("\nEmin_PCG is converged in iter: %d , diffE = %10.7e\n", iter, diffE)
            end
            break
        end

        g_old = copy(g)
        d_old = copy(d)
        Kg_old = copy(Kg)
        Etot_old = Etot

        flush(stdout)
    end

    #GC.gc(true)
    println()
    CuArrays.memory_status()

    # Calculate eigenvalues
    Hr = CuArrays.zeros(ComplexF64, Nstates, Nstates)
    for ispin = 1:Nspin, ik = 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        ortho_gram_schmidt!( psiks[ikspin] )
        Hr[:] = psiks[ikspin]' * op_H(Ham, psiks[ikspin]) # no need to convert to Hermitian
        Ham.electrons.ebands_gpu[:,ik], evecs = eigen_heevd_gpu(Hr)
        psiks[ikspin] = psiks[ikspin]*evecs
        # XXX preallocate evecs?
        #Hr_cpu = collect(Hr)
        #println("Hr_cpu")
        #println(Hr_cpu)
        #evals, evecs = eigen( Hermitian(Hr_cpu) )  # do diagonalization at CPU
        #Ham.electrons.ebands[:,ik] = evals
        #psiks[ikspin] = psiks[ikspin]*CuArray(evecs) # duh!!!
    end

    if verbose && print_final_ebands
        @printf("\n")
        @printf("----------------------------\n")
        @printf("Final Kohn-Sham eigenvalues:\n")
        @printf("----------------------------\n")
        @printf("\n")
        Ham.electrons.ebands[:] = collect( Ham.electrons.ebands_gpu[:] )
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

    return

end
