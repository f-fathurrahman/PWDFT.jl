function linmin_grad_vec!( Ham::Hamiltonian, psiks::BlochWavefunc, g, d, α_t::Float64, α )

    psic = zeros_BlochWavefunc(Ham)
    gt = zeros_BlochWavefunc(Ham)

    Nkspin = length(psiks)
    for i in 1:Nkspin
        psic[i] = psiks[i] + α_t*d[i]
        ortho_sqrt!( psic[i] )
    end

    calc_grad!( Ham, psic, gt )

    for i in 1:Nkspin
        denum = 2.0*real( dot(g[i] - gt[i], d[i]) )
        if denum != 0.0
            α[i] = abs( α_t * 2.0*real( dot(g[i], d[i]) )/denum )
        else
            α[i] = 0.0
        end
    end
    return
end


function KS_solve_Emin_PCG_vec!(
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

    d = deepcopy(Kg)

    # Constrain
    constrain_search_dir!( d, psiks )

    gPrevUsed = true

    α_t = 3e-5
    α = zeros(Float64,Nkspin)
    β = zeros(Float64,Nkspin)
    dotgd = zeros(Float64,Nkspin)
    gKnorm = zeros(Float64,Nkspin)
    gKnormPrev = zeros(Float64,Nkspin)
    dotgPrevKg = zeros(Float64,Nkspin)

    force_grad_dir = true

    Etot_old = Etot
    Nconverges = 0
    
    cg_test = 0.0

    if verbose
        @printf("\n")
        @printf("Minimizing Kohn-Sham energy using PCG vec\n")
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

        for i in 1:Nkspin
            gKnorm[i] = 2*real(dot(g[i], Kg[i]))
        end
        
        if !force_grad_dir
            
            for i in 1:Nkspin
                dotgd[i] = 2*real(dot(g[i], d[i]))
            end


            if gPrevUsed
                for i in 1:Nkspin
                    dotgPrevKg[i] = 2*real(dot(gPrev[i], Kg[i]))
                end
            else
                dotgPrevKg .= 0.0
            end

            for i in 1:Nkspin
                β[i] = (gKnorm[i] - dotgPrevKg[i])/gKnormPrev[i] # Polak-Ribiere
                #β = gKnorm/gKnormPrev # Fletcher-Reeves
                #β = (gKnorm - dotgPrevKg) / ( dotgd - dot_BlochWavefunc(d,gPrev) )
                #β = gKnorm/dot_BlochWavefunc(g .- gPrev, d_old)
                if β[i] < 0.0
                    #println("Resetting β")
                    β[i] = 0.0
                end
            end

        end # force_grad_dir
        #println("β = ", β)

        force_grad_dir = false
        
        if gPrevUsed
            gPrev = deepcopy(g)
        end
        gKnormPrev[:] = gKnorm[:]

        # Update search direction
        for i in 1:Nkspin
            d_old[i] = copy(d[i])
            d[i] = -Kg[i] + β[i]*d[i]
        end

        constrain_search_dir!( d, psiks )

        linmin_grad_vec!( Ham, psiks, g, d, α_t, α )
        #println("α = ", α)

        Rhoe_old = copy(Ham.rhoe)
        # Update psiks, Rhoe, etc
        do_step!( psiks, α, d )
        #println("α = ", α)
        Etot = calc_energies_grad!( Ham, psiks, g, Kg )
        
        diffE = Etot_old - Etot
        #norm_g = sqrt(2.0*real(dot(g,Kg))/length(g))
        norm_g = norm(g)/length(g)
        #norm_g = 2*real(dot(g,g))
        mae_rhoe = sum( abs.( Ham.rhoe - Rhoe_old ) )/(Npoints*Nspin)
        @printf("Emin_PCG_vec step %8d = %18.10f  %12.7e %12.7e %12.7e\n", iter, Etot, diffE, norm_g, mae_rhoe)
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
