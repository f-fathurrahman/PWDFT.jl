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

    Nconverged = 0
    # Begin iter
    for iter in 1:NiterMax

        @printf("---------------------\n")
        @printf("Begin iteration #%4d\n", iter)
        @printf("---------------------\n")

        ss = dot(evars.psiks[1], g.psiks[1])
        @printf("dot evars.psi and g.psi     = %18.10f\n", real(ss))

        ss = dot(diagm(0 => Ham.electrons.ebands[:,1]), g.Haux[1]) # ikspin=1
        @printf("dot diagm(eorbs) and g.Haux = %18.10f\n", real(ss))

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
        if α > 2.0
            println("α is too large: ", α)
            println("Limiting it to 2.0")
            α = 2.0
        end
        if α_Haux > 2.0
            println("α_Haux is too large: ", α_Haux)
            println("Limiting it to 2.0")
            α_Haux= 2.0
        end
        println("α      = ", α)
        println("α_Haux = ", α_Haux)
        do_step!( Ham, α, α_Haux, evars, d, subrot )

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
        if abs(diffE) < 1e-6
            Nconverged = Nconverged + 1
        else
            Nconverged = 0
        end
        if Nconverged >= 2
            println("Emin_SD_Haux is converged")
            break
        end
        calc_Hsub_eigs!(evars)
        print_ebands_Hsub_eigs(Ham, evars)
    end

    println()
    println(Ham.energies)

    return
end