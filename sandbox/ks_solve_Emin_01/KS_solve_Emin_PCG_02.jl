"""
Solves Kohn-Sham problem using direct energy minimization as described
by Ismail-Beigi and Arias.
"""
function KS_solve_Emin_PCG_02!(
    Ham::Hamiltonian;
    startingwfc=:random,
    savewfc=false,
    startingrhoe=:gaussian,
    skip_initial_diag=false,
    α_t=3e-5,
    NiterMax=200,
    verbose=true,
    print_final_ebands=false,
    print_final_energies=true,
    i_cg_beta=2,
    etot_conv_thr=1e-6,
    kT=0.001
)

    pw = Ham.pw
    electrons = Ham.electrons
    
    Nstates = electrons.Nstates
    Nelectrons = electrons.Nelectrons
    
    Ns = pw.Ns
    Npoints = prod(Ns)
    CellVolume = pw.CellVolume
    Ngw = pw.gvecw.Ngw
    Ngwx = pw.gvecw.Ngwx
    Nkpt = pw.gvecw.kpoints.Nkpt
    wk = pw.gvecw.kpoints.wk

    Nspin = electrons.Nspin
    Nkspin = Nkpt*Nspin

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

    Hsub = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    for ikspin = 1:Nkspin
        Hsub[ikspin] = zeros(ComplexF64,Nstates,Nstates)
    end
    Δ_Haux      = copy(Hsub)
    Δt_Haux     = copy(Hsub)
    Δ_Haux_old  = copy(Hsub)
    g_Haux      = copy(Hsub)
    gt_Haux     = copy(Hsub)    
    Kg_Haux     = copy(Hsub)
    g_Haux_old  = copy(Hsub)
    Kg_Haux_old = copy(Hsub)    
    Haux_c      = copy(Hsub)
    d_Haux      = copy(Hsub)
    d_Haux_old  = copy(Hsub)
    U_Haux       = copy(Hsub)

    Haux = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    eta = zeros(Float64,Nstates,Nkspin)
    for ikspin = 1:Nkspin
        Haux[ikspin] = rand(ComplexF64, Nstates, Nstates)
        Haux[ikspin] = 0.5*(Haux[ikspin] + Haux[ikspin]')
        eta[:,ikspin], U_Haux[ikspin] = eigen(Haux[ikspin])
        Haux[ikspin] = U_Haux[ikspin]*Haux[ikspin]*U_Haux[ikspin]'
    end

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

    # Rotate
    for i in 1:Nkspin
        psiks[i] = psiks[i]*U_Haux[i]
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
    
    # Symmetrize Rhoe if needed
    if Ham.sym_info.Nsyms > 1
        symmetrize_rhoe!( Ham, rhoe_symmetrizer, Rhoe )
    end

    update!(Ham, Rhoe)

    evals = zeros(Nstates,Nkspin)



    β_Haux = zeros(Nkspin)
    α_Haux = zeros(Nkspin)

    α_t_Haux = 1e-3

    Etot_old = 0.0

    Ham.electrons.Focc[:,:], E_fermi = calc_Focc( Nelectrons, wk, kT, eta, Nspin )
    Entropy = calc_entropy( wk, kT, eta, E_fermi, Nspin )
    
    println("E_fermi = ", E_fermi)
    for ist = 1:Nstates
        println(eta[ist,1], " ", Ham.electrons.Focc[ist,1])
    end

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( Ham.atoms )

    # calculate PspCore energy
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )

    # Calculate energy at this psi
    energies = calc_energies(Ham, psiks)
    Ham.energies = energies
    Etot = sum(energies)

    Nconverges = 0

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

    Kscalar = 0.05

    for iter = 1:NiterMax

        for ispin = 1:Nspin
        for ik = 1:Nkpt

            Ham.ik = ik
            Ham.ispin = ispin
            ikspin = ik + (ispin - 1)*Nkpt

            g[ikspin], g_Haux[ikspin], _, Δ_Haux[ikspin] =
            calc_grad_Haux( Ham, psiks[ikspin], eta[:,ikspin], kT )
            
            Kg[ikspin] = Kprec( Ham.ik, pw, g[ikspin] )
            
            Kg_Haux[ikspin] = Kscalar*g_Haux[ikspin]

            # XXX: define function trace for real(sum(conj(...)))
            if iter != 1
                
                if i_cg_beta == 1
                    β[ikspin] =
                    real(sum(conj(g[ikspin]).*Kg[ikspin]))/real(sum(conj(g_old[ikspin]).*Kg_old[ikspin]))
                
                elseif i_cg_beta == 2
                    
                    β[ikspin] =
                    real(sum(conj(g[ikspin]-g_old[ikspin]).*Kg[ikspin]))/real(sum(conj(g_old[ikspin]).*Kg_old[ikspin]))
                    
                    #
                    #num1 = real(sum(conj(g_Haux[ikspin]-g_Haux_old[ikspin]).*Kg_Haux[ikspin]))
                    #denum1 = real(sum(conj(g_Haux_old[ikspin]).*Kg_Haux_old[ikspin]))
                    
                    num1 = real(sum( conj(g_Haux[ikspin]) .* Δ_Haux[ikspin] ) )
                    denum1 = real(sum( conj(g_Haux_old[ikspin]) .* Δ_Haux_old[ikspin] ) )
                    β_Haux[ikspin] = num1/denum1

                    println("num1   = ", num1)
                    println("denum1 = ", denum1)
                    println("β Haux = ", β_Haux[ikspin])
                
                elseif i_cg_beta == 3
                    β[ikspin] =
                    real(sum(conj(g[ikspin]-g_old[ikspin]).*Kg[ikspin]))/real(sum(conj(g[ikspin]-g_old[ikspin]).*d[ikspin]))
                
                elseif i_cg_beta == 4
                    β[ikspin] =
                    real(sum(conj(g[ikspin]).*Kg[ikspin]))/real(sum((g[ikspin]-g_old[ikspin]).*conj(d_old[ikspin])))
                
                else
                    β[ikspin] = 0.0
                    β_Haux[ikspin] = 0.0

                end
            end
            
            if β[ikspin] < 0.0
                β[ikspin] = 0.0
            end

            if β_Haux[ikspin] < 0.0
                β_Haux[ikspin] = 0.0
            end


            d[ikspin] = -Kg[ikspin] + β[ikspin] * d_old[ikspin]
            psic[ikspin] = ortho_sqrt(psiks[ikspin] + α_t*d[ikspin])

            d_Haux[ikspin] = -Kg_Haux[ikspin] + β_Haux[ikspin] * d_Haux_old[ikspin]

            Haux_c[ikspin] = Haux[ikspin] + α_t*d_Haux[ikspin]

            Haux_c[ikspin] = 0.5*(Haux_c[ikspin] + Haux_c[ikspin]')
            eta[:,ikspin], U_Haux[ikspin] = eigen(Haux_c[ikspin])
            Haux_c[ikspin] = U_Haux[ikspin]*Haux_c[ikspin]*U_Haux[ikspin]'            

            psic[ikspin] = psic[ikspin]*U_Haux[ikspin] # rotate

        end # ik
        end # ispin


        Ham.electrons.Focc[:,:], E_fermi = calc_Focc( Nelectrons, wk, kT, eta, Nspin )
        Entropy = calc_entropy( wk, kT, eta, E_fermi, Nspin )

        calc_rhoe!( Ham, psic, Rhoe )
        # Symmetrize Rhoe if needed
        if Ham.sym_info.Nsyms > 1
            symmetrize_rhoe!( Ham, rhoe_symmetrizer, Rhoe )
        end
        update!(Ham, Rhoe)

        for ispin = 1:Nspin
        for ik = 1:Nkpt
            
            Ham.ik = ik
            Ham.ispin = ispin
            ikspin = ik + (ispin - 1)*Nkpt

            gt[ikspin], gt_Haux[ikspin], _, Δt_Haux[ikspin] =
            calc_grad_Haux( Ham, psic[ikspin], eta[:,ikspin], kT )

            denum = real(sum(conj(g[ikspin]-gt[ikspin]).*d[ikspin]))
            #println("denum psi = ", denum)
            if denum != 0.0
                α[ikspin] = abs( α_t*real(sum(conj(g[ikspin]).*d[ikspin]))/denum )
            else
                α[ikspin] = 0.0
            end
            # Update wavefunction
            psiks[ikspin] = psiks[ikspin] + α[ikspin]*d[ikspin]
            # and orthonormalize it
            psiks[ikspin] = ortho_sqrt(psiks[ikspin])

            # Update Haux
            #denum = real(sum(conj(g_Haux[ikspin]-gt_Haux[ikspin]).*d_Haux[ikspin]))
            denum = real(sum(conj(Δ_Haux[ikspin]-Δt_Haux[ikspin]).*d_Haux[ikspin]))

            if abs(denum) > 1e-15
                α_Haux[ikspin] = abs( α_t*real(sum(conj(Δ_Haux[ikspin]).*d_Haux[ikspin]))/denum )
            else
                α_Haux[ikspin] = 0.0
                #α_Haux[ikspin] = α_t
            end
            Haux[ikspin] = Haux[ikspin] + α_Haux[ikspin]*d_Haux[ikspin]
            Haux[ikspin] = 0.5*(Haux[ikspin] + Haux[ikspin]')
            eta[:,ikspin], U_Haux[ikspin] = eigen(Haux[ikspin])
            Haux[ikspin] = U_Haux[ikspin]*Haux[ikspin]*U_Haux[ikspin]'            

            psiks[ikspin] = psiks[ikspin]*U_Haux[ikspin] # rotate
            # also rotate Haux ?
        end
        end

        Ham.electrons.Focc[:,:], E_fermi = calc_Focc( Nelectrons, wk, kT, eta, Nspin )
        Entropy = calc_entropy( wk, kT, eta, E_fermi, Nspin )

        # Calculate Hsub (for comparison with Haux)
        for ispin = 1:Nspin
        for ik = 1:Nkpt
            Ham.ik = ik
            Ham.ispin = ispin
            ikspin = ik + (ispin - 1)*Nkpt
            #
            Hsub[ikspin] = psiks[ikspin]' * ( Ham*psiks[ikspin] )

            evalsHsub = eigvals(Hsub[ikspin])

            println("diff Haux = ", sum(Hsub[ikspin] - Haux[ikspin]))
            println("diff evals = ", sum(evalsHsub - eta[:,ikspin]))

        end
        end



        println("\nOccupation numbers")
        println("E_fermi = ", E_fermi)
        for i = 1:Nstates
            @printf("eta, Focc = %18.10f %18.10f\n", eta[i,1], Ham.electrons.Focc[i,1])
        end

        calc_rhoe!( Ham, psiks, Rhoe )
        # Symmetrize Rhoe if needed
        if Ham.sym_info.Nsyms > 1
            symmetrize_rhoe!( Ham, rhoe_symmetrizer, Rhoe )
        end
        update!(Ham, Rhoe)

        Ham.energies = calc_energies( Ham, psiks )
        Ham.energies.mTS = Entropy

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

        g_Haux_old = copy(g_Haux)
        Δ_Haux_old = copy(g_Haux)
        d_Haux_old = copy(d_Haux)
        Kg_Haux_old = copy(Kg_Haux)

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
