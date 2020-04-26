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

    rotPrev = Vector{Matrix{ComplexF64}}(undef,Nkspin)
    rotPrevC = Vector{Matrix{ComplexF64}}(undef,Nkspin)
    rotPrevCinv = Vector{Matrix{ComplexF64}}(undef,Nkspin)
    for i in 1:Nkspin
        rotPrev[i] = diagm( 0 => ones(ComplexF64,Nstates) )
        rotPrevC[i] = diagm( 0 => ones(ComplexF64,Nstates) )
        rotPrevCinv[i] = diagm( 0 => ones(ComplexF64,Nstates) )
    end

    Ham.energies.NN = calc_E_NN(Ham.atoms)
    
    Etot = compute!( Ham, evars, g, Kg, kT, rotPrevCinv, rotPrev )
    Etot_old = Etot

    @printf("Initial energies = %18.10f\n", Etot)

    d = deepcopy(g)

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
            d.psiks[i] = -g.psiks[i]
            d.Haux[i]  = -g.Haux[i] 
        end

        constrain_search_dir!( d, evars )

        α = linmin_grad_v1!( Ham, evars, g, d, kT, rotPrev, rotPrevC, rotPrevCinv )
        println("α     = ", α)

        do_step!( α, evars, d, rotPrev, rotPrevC, rotPrevCinv )

        Etot_old = Etot
        Etot = compute!( Ham, evars, g, Kg, kT, rotPrevCinv, rotPrev )
        #println(Ham.energies)
        diffE = Etot - Etot_old
        @printf("Emin_SD_Haux: %5d %18.10f %18.10e ", iter, Etot, abs(diffE))
        if diffE > 0
            println("Energy is not reducing")
        else
            println()
        end

    end

    println(Ham.energies)

    return
end