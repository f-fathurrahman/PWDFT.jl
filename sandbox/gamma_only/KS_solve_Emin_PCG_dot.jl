function KS_solve_Emin_PCG_dot!(
    Ham::Hamiltonian, psiks::BlochWavefunc;
    startingrhoe=:gaussian,
    skip_initial_diag=false,
    α_t=3e-5, NiterMax=200, verbose=true,
    print_final_ebands=false, print_final_energies=true,
    etot_conv_thr=1e-6
)

    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = length(psiks)

    g = zeros_BlochWavefunc(Ham)
    Kg = zeros_BlochWavefunc(Ham)
    gPrev = zeros_BlochWavefunc(Ham)

    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Nstates = Ham.electrons.Nstates

    setup_guess_wavefunc!( Ham, psiks, startingrhoe, skip_initial_diag )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( Ham.atoms )

    # No need to orthonormalize
    Etot = calc_energies_grad!( Ham, psiks, g, Kg )

    d = deepcopy(Kg)
    # Constrain
    #constrain_search_dir!( d, psiks )

    gPrevUsed = true

    α = 0.0
    β = 0.0
    gKnorm = 0.0
    gKnormPrev = 0.0
    force_grad_dir = true

    Etot_old = Etot
    Nconverges = 0
    
    cg_test = 0.0

    if verbose
        @printf("\n")
        @printf("Minimizing Kohn-Sham energy using PCG dot\n")
        @printf("-----------------------------------------\n")
        @printf("NiterMax  = %d\n", NiterMax)
        @printf("α_t       = %e\n", α_t)
        @printf("conv_thr  = %e\n", etot_conv_thr)
        @printf("\n")
    end

    Rhoe_old = zeros(Float64,Npoints,2)
    kpoints = Ham.pw.gvecw.kpoints

    for iter in 1:NiterMax

        gKnorm = dot_BlochWavefunc(g, Kg)
        
        if !force_grad_dir    
            dotgd = dot_BlochWavefunc(g, d)
            if gPrevUsed
                dotgPrevKg = dot_BlochWavefunc(gPrev, Kg)
            else
                dotgPrevKg = 0.0
            end
            β = (gKnorm - dotgPrevKg)/gKnormPrev # Polak-Ribiere
        end

        if β < 0.0
            println("Resetting β")
            β = 0.0
        end

        force_grad_dir = false
        
        if gPrevUsed
            gPrev = deepcopy(g)
        end
        gKnormPrev = gKnorm

        # Update search direction
        for i in 1:Nkspin
            d[i] = -Kg[i] + β*d[i]
        end

        #constrain_search_dir!( d, psiks )

        α = linmin_grad!( Ham, psiks, g, d, Etot )
        # Limite the value of α if it is too big.
        # At least found in the case of NH3
        if α > 2.0
            α = 2.0
        end

        Rhoe_old = copy(Ham.rhoe)

        # Update psiks
        #do_step!( psiks, α, d )
        for i in 1:Nspin
            psiks[i] = psiks[i] + α*d[i]
            ortho_sqrt!( psiks[i] )
        end
        #println("α = ", α)
        
        # Calculate rhoe, update rhoe, calc energies and grad
        Etot = calc_energies_grad!( Ham, psiks, g, Kg )
        
        diffE = Etot_old - Etot
        #norm_g = 2*real(dot(g,g))
        #mae_rhoe = sum( abs.( Ham.rhoe - Rhoe_old ) )/(Npoints*Nspin)
        @printf("Emin_PCG_dot step %8d = %18.10f  %12.7e\n", iter, Etot, diffE)
        if diffE < 0.0
            println("*** WARNING: Etot is not decreasing")
        end

        if abs(diffE) < etot_conv_thr
            Nconverges = Nconverges + 1
        else
            Nconverges = 0
        end

        if Nconverges >= 2
            @printf("\nEmin_PCG_dot is converged in iter: %d\n", iter)
            break
        end

        Etot_old = Etot

    end

    # Calculate eigenvalues
    evecs = zeros(ComplexF64,Nstates,Nstates)
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        psiks[ikspin] = ortho_gram_schmidt(psiks[ikspin])
        Hr = Hermitian(psiks[ikspin]' * op_H(Ham, psiks[ikspin]))
        evals, evecs = eigen(Hr)
        Ham.electrons.ebands[:,ik] = evals
        psiks[ikspin] = psiks[ikspin]*evecs
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

    return

end
