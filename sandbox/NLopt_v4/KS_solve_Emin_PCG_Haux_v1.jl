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

    rotPrev = Vector{Matrix{ComplexF64}}(undef,Nkspin)
    rotPrevC = Vector{Matrix{ComplexF64}}(undef,Nkspin)
    rotPrevCinv = Vector{Matrix{ComplexF64}}(undef,Nkspin)
    for i in 1:Nkspin
        rotPrev[i] = diagm( 0 => ones(ComplexF64,Nstates) )
        rotPrevC[i] = diagm( 0 => ones(ComplexF64,Nstates) )
        rotPrevCinv[i] = diagm( 0 => ones(ComplexF64,Nstates) )
    end

    Ham.energies.NN = calc_E_NN(Ham.atoms)
    
    #Etot = calc_energies_grad!( Ham, evars, g, Kg, kT )
    Etot = compute!( Ham, evars, g, Kg, kT, rotPrevCinv, rotPrev )
    Etot_old = Etot

    #println(Ham.energies)
    @printf("Initial energies = %18.10f\n", Etot)

    #d = deepcopy(Kg)
    d = deepcopy(g)

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
            #d.psiks[i] = -Kg.psiks[i] + β*d.psiks[i] #-g.psiks[i]
            #d.Haux[i]  = -Kg.Haux[i] + β*d.Haux[i]   #-g.Haux[i] 
            d.psiks[i] = -Kg.psiks[i] + β*d.psiks[i]
            d.Haux[i]  = -Kg.Haux[i] + β_aux*d.Haux[i]
        end

        constrain_search_dir!( d, evars )

        #println("rotPrevCinv")
        #print_vec_mat(rotPrevCinv[1:1])

        #α = linmin_grad!( Ham, evars, g, d, kT, rotPrev, rotPrevC, rotPrevCinv )
        #α = linmin_grad!( Ham, evars.psiks, g.psiks, d.psiks )
        α, α_aux = linmin_grad!( Ham, evars, g, d, kT, rotPrev, rotPrevC, rotPrevCinv )
        #println("α     = ", α)
        #println("α_aux = ", α_aux)
        #α = 3e-5

        #do_step!( α, 0.0, evars, d, rotPrev, rotPrevC, rotPrevCinv )
        #do_step!( 0.0, α, evars, d, rotPrev, rotPrevC, rotPrevCinv )
        #do_step!( α, evars, d, rotPrev, rotPrevC, rotPrevCinv )
        do_step!( α, α_aux, evars, d, rotPrev, rotPrevC, rotPrevCinv )

        Etot_old = Etot
        Etot = compute!( Ham, evars, g, Kg, kT, rotPrevCinv, rotPrev )
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