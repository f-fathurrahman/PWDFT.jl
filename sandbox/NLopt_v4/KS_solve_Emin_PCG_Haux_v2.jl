function KS_solve_Emin_PCG_Haux_v2!( Ham::Hamiltonian, evars::ElecVars;
    etot_conv_thr=1e-6, skip_initial_diag=false, startingrhoe=:gaussian, NiterMax=5, kT=0.01
)

    Nstates = Ham.electrons.Nstates
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin

    g = ElecGradient(Ham)
    Kg = ElecGradient(Ham)
    gPrev = ElecGradient(Ham)

    subrot = SubspaceRotations(Nkspin, Nstates)

    Ham.energies.NN = calc_E_NN(Ham.atoms)
    
    Etot = compute!( Ham, evars, g, Kg, kT, subrot )
    Etot_old = Etot

    @printf("Initial energies = %18.10f\n", Etot)

    d = deepcopy(Kg)

    # Constrain
    constrain_search_dir!( d, evars )

    β = 0.0
    β_aux = 0.0

    gPrevUsed = true
    gKnormPrev = 0.0
    gKnormPrev_aux = 0.0

    force_grad_dir = true

    # Begin iter
    for iter in 1:NiterMax

        #gKnorm = dot_ElecGradient(g, Kg)
        gKnorm, gKnorm_aux = dot_ElecGradient_v2(g, Kg)
        @printf("gKnorm = %18.10e, gKnorm_aux = %18.10e\n", gKnorm, gKnorm_aux)

        if !force_grad_dir
            
            dotgd, dotgd_aux = dot_ElecGradient_v2(g, d)
            
            if gPrevUsed
                dotgPrevKg, dotgPrevKg_aux = dot_ElecGradient_v2(gPrev, Kg)
            else
                dotgPrevKg = 0.0
                dotgPrevKg_aux = 0.0
            end
            β = (gKnorm - dotgPrevKg)/gKnormPrev # Polak-Ribiere
            println("β = ", β)
            if β < 0.0
                println("Resetting β")
                β = 0.0
            end
            β_aux = (gKnorm_aux - dotgPrevKg_aux)/gKnormPrev_aux # Polak-Ribiere
            println("β_aux = ", β_aux)
            if β_aux < 0.0
                println("Resetting β_aux")
                β_aux = 0.0
            end

        end

        force_grad_dir = false

        # Check convergence here ....

        # No convergence yet, continuing ...
        
        if gPrevUsed
            gPrev = deepcopy(g)
        end
        gKnormPrev = copy(gKnorm) # Need copy here?
        gKnormPrev_aux = copy(gKnorm_aux) # Need copy here?

        # Update search direction
        for i in 1:Nkspin
            d.psiks[i] = -Kg.psiks[i] + β*d.psiks[i]
            d.Haux[i]  = -Kg.Haux[i] + β_aux*d.Haux[i]
        end

        constrain_search_dir!( d, evars )

        α, α_aux = linmin_grad_v2!( Ham, evars, g, d, kT, subrot )

        do_step!( Ham, α, α_aux, evars, d, subrot )

        Etot_old = Etot
        Etot = compute!( Ham, evars, g, Kg, kT, subrot )
        #println(Ham.energies)
        diffE = Etot - Etot_old
        println()
        @printf("Emin_PCG_Haux_v2: %5d %18.10f %18.10e ", iter, Etot, abs(diffE))
        if diffE > 0
            println("Energy is not reducing")
        else
            println()
        end
        print_ebands(Ham.electrons, Ham.pw.gvecw.kpoints)


    end

    println(Ham.energies)

    return
end