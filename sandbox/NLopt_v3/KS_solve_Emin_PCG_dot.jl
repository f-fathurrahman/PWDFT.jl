function KS_solve_Emin_PCG_dot!(
    Ham::Hamiltonian, psiks::BlochWavefunc;
    startingrhoe=:gaussian,
    skip_initial_diag=false,
    α_t=3e-5, NiterMax=200, verbose=true,
    print_final_ebands=false, print_final_energies=true,
    i_cg_beta=2, etot_conv_thr=1e-6
)

    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = length(psiks)

    g = zeros_BlochWavefunc(Ham)
    Kg = zeros_BlochWavefunc(Ham)
    gPrev = zeros_BlochWavefunc(Ham)
    d_old = zeros_BlochWavefunc(Ham) # needed for β Dai-Yuan

    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Nstates = Ham.electrons.Nstates

    setup_guess_wavefunc!( Ham, psiks, startingrhoe, skip_initial_diag )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( Ham.atoms )

    # No need to orthonormalize
    Etot = calc_energies_grad!( Ham, psiks, g, Kg )
    #println("Initial Etot = ", Etot)
    #println("Initial dot_BlochWavefunc(g,g) = ", dot_BlochWavefunc(g,g))

    d = deepcopy(Kg)

    # Constrain
    constrain_search_dir!( d, psiks )

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

    Rhoe_old = zeros(Float64,Npoints,2)
    kpoints = Ham.pw.gvecw.kpoints

    for iter in 1:NiterMax

        gKnorm = dot_BlochWavefunc(g, Kg)
        #gKnorm = dot_BlochWavefunc(kpoints, g, Kg)
        
        if !force_grad_dir
            
            dotgd = dot_BlochWavefunc(g, d)
            if gPrevUsed
                dotgPrevKg = dot_BlochWavefunc(gPrev, Kg)
                #dotgPrevKg = dot_BlochWavefunc(kpoints, gPrev, Kg)
            else
                dotgPrevKg = 0.0
            end

            β = (gKnorm - dotgPrevKg)/gKnormPrev # Polak-Ribiere
            #β = gKnorm/gKnormPrev # Fletcher-Reeves
            #β = (gKnorm - dotgPrevKg) / ( dotgd - dot_BlochWavefunc(d,gPrev) )
            #β = gKnorm/dot_BlochWavefunc(g .- gPrev, d_old)
            #β = 0.0
            #denum = sqrt( dot_BlochWavefunc(g,g) * dot_BlochWavefunc(d,d) )
            #println("linmin test: ", dotgd/denum )
            #if gPrevUsed
            #    cg_test  = dotgPrevKg/sqrt(gKnorm*gKnormPrev)
            #    println("CG test: ", cg_test)
            #end
        end

        if β < 0.0
            println("Resetting β")
            β = 0.0
        end

        #println("β = ", β)

        force_grad_dir = false
        
        if gPrevUsed
            gPrev = deepcopy(g)
        end
        gKnormPrev = gKnorm

        # Update search direction
        for i in 1:Nkspin
            d_old[i] = copy(d[i])
            d[i] = -Kg[i] + β*d[i]
        end

        constrain_search_dir!( d, psiks )

        _, α = linmin_grad!( Ham, psiks, g, d, Etot )
        #println("α = ", α)

        Rhoe_old = copy(Ham.rhoe)
        # Update psiks
        do_step!( psiks, α, d )
        #println("α = ", α)
        # Calculate rhoe, update rhoe, calc energies and grad
        Etot = calc_energies_grad!( Ham, psiks, g, Kg )
        
        diffE = Etot_old - Etot
        #norm_g = sqrt(2.0*real(dot(g,Kg))/length(g))
        norm_g = norm(g)/length(g)
        #norm_g = 2*real(dot(g,g))
        mae_rhoe = sum( abs.( Ham.rhoe - Rhoe_old ) )/(Npoints*Nspin)
        @printf("Emin_PCG_dot step %8d = %18.10f  %12.7e %12.7e %12.7e\n", iter, Etot, diffE, norm_g, mae_rhoe)
        if diffE < 0.0
            println("*** WARNING: Etot is not decreasing")
        end

        if abs(diffE) < etot_conv_thr
            Nconverges = Nconverges + 1
        else
            Nconverges = 0
        end

        #if (Nconverges >= 2) && (norm_g > etot_conv_thr*5) # (mae_rhoe >= etot_conv_thr) #(norm_g >= etot_conv_thr/10)
        #    println("Probably early convergence, continuing ...")
        #    Nconverges = 0
        #end
        
        #if (Nconverges >= 2) && (2*real(dot(g,g)) > 1e-10) #(norm_g > 1e-6)
        #    println("Probably early convergence, continuing ...")
        #    Nconverges = 0
        #end

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
        psiks[ikspin] = ortho_sqrt(psiks[ikspin])
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
