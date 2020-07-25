function KS_solve_Emin_PCG_dot!(
    Ham::HamiltonianGamma, psis::BlochWavefuncGamma;
    startingrhoe=:random,
    skip_initial_diag=true,
    α_t=3e-5, NiterMax=200, verbose=true,
    print_final_ebands=false, print_final_energies=true,
    i_cg_beta=2, etot_conv_thr=1e-6
)

    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Nstates = Ham.electrons.Nstates

    Rhoe = zeros(Float64,Npoints,Nspin)
    Rhoe_old = zeros(Float64,Npoints,Nspin)

    g     = zeros_BlochWavefuncGamma(Ham)
    Kg    = zeros_BlochWavefuncGamma(Ham)
    gPrev = zeros_BlochWavefuncGamma(Ham)

    Hsub = Vector{Matrix{ComplexF64}}(undef,Nspin)
    for ispin in 1:Nspin
        Hsub[ispin] = zeros(ComplexF64,Nstates,Nstates)
    end

    setup_guess_wavefunc!( Ham, psis, startingrhoe, skip_initial_diag )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( Ham.atoms )

    println()
    println("Initial dot(psis,psis) = ", 2*real(dot_BlochWavefuncGamma(psis,psis)))

    # No need to orthonormalize
    Etot = calc_energies_grad!( Ham, psis, g, Kg, Hsub )

    println("Initial dot(psis,psis) = ", 2*real(dot_BlochWavefuncGamma(psis,psis)))
    
    println("Initial Etot = ", Etot)
    println("Initial dot(g,g) = ", 2*real(dot_BlochWavefuncGamma(g,g)))
    println("Initial dot(Kg,Kg) = ", 2*real(dot_BlochWavefuncGamma(Kg,Kg)))

    d = deepcopy(Kg)
    
    println("Before constrain_search_dir dot(d,d) = ", 2*real(dot_BlochWavefuncGamma(d,d)))

    # Constrain
    constrain_search_dir!( d, psis )

    println("After constrain_search_dir dot(d,d) = ", 2*real(dot_BlochWavefuncGamma(d,d)))

    println()

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

    for iter in 1:NiterMax

        gKnorm = 2*real(dot_BlochWavefuncGamma(g, Kg))
        
        if !force_grad_dir    
            dotgd = 2*real(dot_BlochWavefuncGamma(g, d))
            if gPrevUsed
                dotgPrevKg = 2*real(dot_BlochWavefuncGamma(gPrev, Kg))
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
        for i in 1:Nspin
            d.data[i] = -Kg.data[i] + β*d.data[i]
        end

        constrain_search_dir!( d, psis )

        α = linmin_grad!( Ham, psis, g, d )
        #println("α = ", α)

        Rhoe_old = copy(Ham.rhoe)
        
        # Update psis
        for i in 1:Nspin
            psis.data[i] = psis.data[i] + α*d.data[i]
            ortho_GS_gamma!( psis.data[i] )
        end

        Etot = calc_energies_grad!( Ham, psis, g, Kg, Hsub )
        
        diffE = Etot_old - Etot
        #norm_g = sqrt(2.0*real(dot(g,Kg))/length(g))
        #norm_g = norm(g.data)/length(g.data)
        #norm_g = 2*real(dot_BlochWavefuncGamma(g,g))
        #mae_rhoe = sum( abs.( Ham.rhoe - Rhoe_old ) )/(Npoints*Nspin)
        #@printf("Emin_PCG_dot gamma step %8d = %18.10f  %12.7e %12.7e %12.7e\n", iter, Etot, diffE, norm_g, mae_rhoe)
        @printf("Emin_PCG_dot gamma step %8d = %18.10f  %12.7e\n", iter, Etot, diffE)
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
    for ispin in 1:Nspin
        evals, evecs = eigen(Hermitian(Hsub[ispin]))
        Ham.electrons.ebands[:,ispin] = evals
        psis.data[ispin] = psis.data[ispin]*evecs
    end

    #if verbose && print_final_ebands
    #    @printf("\n")
    #    @printf("----------------------------\n")
    #    @printf("Final Kohn-Sham eigenvalues:\n")
    #    @printf("----------------------------\n")
    #    @printf("\n")
    #    print_ebands(Ham.electrons, Ham.pw.gvecw.kpoints, unit="eV")
    #end

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
