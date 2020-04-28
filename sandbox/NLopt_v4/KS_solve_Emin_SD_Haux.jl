function KS_solve_Emin_SD_Haux!( Ham::Hamiltonian, evars::ElecVars;
    etot_conv_thr=1e-6, skip_initial_diag=false, startingrhoe=:gaussian, NiterMax=100, kT=0.01
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

    d = deepcopy(g) # XXX should only allocate memory

    # Constrain
    constrain_search_dir!( d, evars )

    β = 0.0

    # Begin iter
    for iter in 1:NiterMax

        gKnorm = dot_ElecGradient(g, Kg)
        @printf("gKnorm = %18.10e\n", gKnorm)

        # Check convergence here ....

        # No convergence yet, continuing ...

        # Update search direction
        for i in 1:Nkspin
            d.psiks[i] = -Kg.psiks[i]
            d.Haux[i]  = -Kg.Haux[i] 
        end

        constrain_search_dir!( d, evars )

        α, α_Haux = linmin_grad_v2!( Ham, evars, g, d, kT, subrot )
        println("α      = ", α)
        println("α_Haux = ", α_Haux)
        do_step!( α, α_Haux, evars, d, subrot )

        Etot_old = Etot
        Etot = compute!( Ham, evars, g, Kg, kT, subrot )
        #println(Ham.energies)
        diffE = Etot - Etot_old
        @printf("Emin_SD_Haux: %5d %18.10f %18.10e ", iter, Etot, abs(diffE))
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