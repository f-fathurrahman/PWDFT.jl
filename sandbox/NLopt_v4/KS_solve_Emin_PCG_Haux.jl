function KS_solve_Emin_PCG_Haux!(
     Ham::Hamiltonian, evars::ElecVars;
    etot_conv_thr=1e-6, skip_initial_diag=false, startingrhoe=:gaussian, NiterMax=5, kT=0.01
)

    g = ElecGradient(Ham)
    Kg = ElecGradient(Ham)
    gPrev = ElecGradient(Ham)

    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nstates = Ham.electrons.Nstates
    Nkspin = Nkpt*Nspin

    subrot = SubspaceRotations(Nkspin, Nstates)

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( Ham.atoms )

    Etot = compute!( Ham, evars, g, Kg, kT, subrot )

    d = deepcopy(Kg)

    # Constrain
    constrain_search_dir!( d, evars )

    gPrevUsed = true

    minim_params = MinimizeParams()

    αt_start = minim_params.αt_start
    αt_min = minim_params.αt_min
    updateTestStepSize = minim_params.updateTestStepSize
    αt = αt_start
    α = αt

    β = 0.0
    β_aux = 0.0
    
    gKnorm = 0.0
    gKnorm_aux = 0.0
    
    gKnormPrev = 0.0
    gKnormPrev_aux = 0.0

    force_grad_dir = true

    Etot_old = Etot
    Nconverges = 0
    
    cg_test = 0.0
    norm_g = 0.0

    for iter in 1:NiterMax

        println("\nBegin iter = ", iter)

        #gKnorm = dot_ElecGradient(g,Kg)
        #if !force_grad_dir
        #    dotgd = dot_ElecGradient(g, d)
        #    if gPrevUsed
        #        dotgPrevKg = dot_ElecGradient(gPrev, Kg)
        #    else
        #        dotgPrevKg = 0.0
        #    end
        #    β = (gKnorm - dotgPrevKg)/gKnormPrev # Polak-Ribiere
        #    if β < 0.0
        #        println("Resetting β")
        #        β = 0.0
        #    end
        #end

        #gKnorm = dot_ElecGradient(g, Kg)
        
        #gKnorm, gKnorm_aux = dot_ElecGradient_v2(g, Kg)
        #@printf("gKnorm = %18.10e, gKnorm_aux = %18.10e\n", gKnorm, gKnorm_aux)

        #if !force_grad_dir
        #    
        #    dotgd, dotgd_aux = dot_ElecGradient_v2(g, d)
        #    
        #    if gPrevUsed
        #        dotgPrevKg, dotgPrevKg_aux = dot_ElecGradient_v2(gPrev, Kg)
        #    else
        #        dotgPrevKg = 0.0
        #        dotgPrevKg_aux = 0.0
        #    end
        #    β = (gKnorm - dotgPrevKg)/gKnormPrev # Polak-Ribiere
        #    println("β = ", β)
        #    if β < 0.0
        #        println("Resetting β")
        #        β = 0.0
        #    end
        #    β_aux = (gKnorm_aux - dotgPrevKg_aux)/gKnormPrev_aux # Polak-Ribiere
        #    println("β_aux = ", β_aux)
        #    if β_aux < 0.0
        #        println("Resetting β_aux")
        #        β_aux = 0.0
        #    end
        #end


        force_grad_dir = false

        # Check convergence here ....

        # No convergence yet, continuing ...
        
        if gPrevUsed
            gPrev = deepcopy(g)
        end
        gKnormPrev = copy(gKnorm)
        gKnormPrev_aux = copy(gKnorm_aux)

        # Update search direction
        for i in 1:Nkspin
            d.psiks[i] = -Kg.psiks[i] + β*d.psiks[i]
            d.Haux[i]  = -Kg.Haux[i] + β*d.Haux[i]
            #d.psiks[i] = -Kg.psiks[i] + β*d.psiks[i]
            #d.Haux[i]  = -Kg.Haux[i] + β_aux*d.Haux[i]            
        end

        constrain_search_dir!( d, evars )

        # Line minimization
        linmin_success, α, αt = linmin_quad!( Ham, evars, g, d, kT, subrot, α, αt, Etot, minim_params )

        println("linmin_success = ", linmin_success)
        @printf("α = %18.10e, αt = %18.10e\n", α, αt)

        #linmin_success, α = linmin_armijo!( Ham, psiks, g, d, Etot )

        # Using alternative line minimization
        #linmin_success, α = linmin_grad!( Ham, psiks, g, d, Etot )

        if linmin_success
            #
            do_step!( Ham, α, evars, d, subrot )
            Etot = compute!( Ham, evars, g, Kg, kT, subrot )
            #
            if updateTestStepSize
                αt = α
                if αt < αt_min    # bad step size
                    @printf("Bad step size: αt is reset to αt_start: %f\n", αt_start)
                    αt = αt_start # make sure next test step size is not too bad
                end
            end
        else
            # linmin failed:
            do_step!( Ham, -α, evars, d, subrot )
            Etot = compute!( Ham, evars, g, Kg, kT, subrot )
            
            @printf("linmin is failed: Update psiks by αt_min = %e, Etot = %18.10f\n", αt_min, Etot)

            if β >= 1e-10   # should be compared with small number
                # Failed, but not along the gradient direction:
                @printf("Step failed: resetting search direction.\n")
                forceGradDirection = true # reset search direction
            else
                # Failed along the gradient direction
                @printf("Step failed along negative gradient direction.\n")
                @printf("Probably at roundoff error limit. (Should stop)\n")
                #return
            end
        end
        
        #norm_g = norm(g) #/(Nkspin*Nstates)
        #norm_g = 2*real(dot(g,g))
        diffE = Etot_old - Etot
        @printf("Emin_PCG_Haux step %8d = %18.10f  %10.7e\n", iter, Etot, diffE)
        if diffE < 0.0
            println("*** WARNING: Etot is not decreasing")
        end
        calc_Hsub_eigs!( evars )
        print_ebands_Hsub_eigs(Ham, evars)

        #if abs(diffE) < etot_conv_thr
        #    Nconverges = Nconverges + 1
        #else
        #    Nconverges = 0
        #end
        #
        #if Nconverges >= 2
        #    @printf("\nEmin_PCG is converged in iter: %d\n", iter)
        #    break
        #end

        Etot_old = Etot

    end

    println(Ham.energies)

    println("Leaving KS_solve_Emin_PCG_Haux")

end


