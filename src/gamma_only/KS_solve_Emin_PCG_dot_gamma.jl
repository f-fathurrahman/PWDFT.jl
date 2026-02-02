function KS_solve_Emin_PCG_dot!(
    Ham::HamiltonianGamma, psis::BlochWavefuncGamma;
    startingrhoe=:random,
    skip_initial_diag=true,
    α_t=3e-5, NiterMax=200, verbose=true,
    print_final_ebands=false, print_final_energies=true,
    etot_conv_thr=1e-6
)
    # Only Polak-Ribiere form is implemented

    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin_wf
    Nstates = Ham.electrons.Nstates

    Rhoe = zeros(Float64,Npoints,Nspin)
    Rhoe_old = zeros(Float64,Npoints,Nspin)

    g     = zeros_BlochWavefuncGamma(Ham)
    Kg    = zeros_BlochWavefuncGamma(Ham)
    gPrev = zeros_BlochWavefuncGamma(Ham)
    psic = zeros_BlochWavefuncGamma(Ham)
    gt = zeros_BlochWavefuncGamma(Ham)

    Hsub = Vector{Matrix{ComplexF64}}(undef,Nspin)
    for ispin in 1:Nspin
        Hsub[ispin] = zeros(ComplexF64,Nstates,Nstates)
    end

    setup_guess_wavefunc!( Ham, psis, startingrhoe, skip_initial_diag )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( Ham.atoms )

    # No need to orthonormalize
    Etot = calc_energies_grad!( Ham, psis, g, Kg, Hsub )

    d = deepcopy(Kg)
    # Constrain
    constrain_search_dir!( d, psis )

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

    for iter in 1:NiterMax

        gKnorm = real(dot(g, Kg))
        
        if !force_grad_dir    
            dotgd = real(dot(g, d))
            if gPrevUsed
                dotgPrevKg = real(dot(gPrev, Kg))
            else
                dotgPrevKg = 0.0
            end
            β = (gKnorm - dotgPrevKg)/gKnormPrev # Polak-Ribiere
        end

        if β < 0.0
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
        α = linmin_grad!( Ham, psis, g, d, psic, gt, αt=1e-5 )
        
        # Update psis
        for i in 1:Nspin
            psis.data[i] = psis.data[i] + α*d.data[i]
            ortho_sqrt_gamma!( psis.data[i] )
        end

        Etot = calc_energies_grad!( Ham, psis, g, Kg, Hsub )
        
        diffE = Etot_old - Etot
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
