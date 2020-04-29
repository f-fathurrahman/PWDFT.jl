function KS_solve_Emin_PCG_Haux_v1!( Ham::Hamiltonian, evars::ElecVars;
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

    gPrevUsed = true
    gKnormPrev = 0.0

    force_grad_dir = true

    # Begin iter
    for iter in 1:NiterMax

        gKnorm = dot_ElecGradient(g, Kg)
        @printf("gKnorm = %18.10e\n", gKnorm)

        if !force_grad_dir
            
            dotgd = dot_ElecGradient(g, d)
            
            if gPrevUsed
                dotgPrevKg = dot_ElecGradient(gPrev, Kg)
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
        end

        force_grad_dir = false

        # Check convergence here ....

        # No convergence yet, continuing ...
        
        if gPrevUsed
            gPrev = deepcopy(g)
        end
        gKnormPrev = copy(gKnorm) # Need copy here?

        # Update search direction
        for i in 1:Nkspin
            d.psiks[i] = -Kg.psiks[i] + β*d.psiks[i] #-g.psiks[i]
            d.Haux[i]  = -Kg.Haux[i] + β*d.Haux[i]   #-g.Haux[i] 
        end

        constrain_search_dir!( d, evars )

        α = linmin_grad_v1!( Ham, evars, g, d, kT, subrot )
        println("α     = ", α)

        do_step!( Ham, α, evars, d, subrot )

        Etot_old = Etot
        Etot = compute!( Ham, evars, g, Kg, kT, subrot )
        #println(Ham.energies)
        diffE = Etot - Etot_old
        @printf("Emin_PCG: %5d %18.10f %18.10e ", iter, Etot, abs(diffE))
        if diffE > 0
            println("Energy is not reducing")
        else
            println()
        end

    end

    println(Ham.energies)

    return
end