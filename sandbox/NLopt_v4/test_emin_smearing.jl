using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("smearing.jl")
include("create_Ham.jl")
include("ElecVars.jl")
include("test_ElecVars.jl")
include("calc_energies_grad.jl")
include("emin_smearing.jl")
include("linmin_grad.jl")

function print_vec_mat( v::Vector{Matrix{ComplexF64}} )
    Nkspin = length(v)
    for i in 1:Nkspin
        println("Real part of ikspin = ", i)
        display(real(v[i])); println()
        println("Imag part of ikspin = ", i)
        display(imag(v[i])); println()    
    end
end

function main()

    Random.seed!(1234)

    kT = 0.01
    #Ham = create_Ham_atom_Al_smearing()
    Ham = create_Ham_Al_fcc_smearing()
    #Ham = create_Ham_atom_Pt_smearing()
    #Ham = create_Ham_Pt_fcc_smearing()
    #println(Ham)

    #test_ElecVars(Ham)

    Nstates = Ham.electrons.Nstates
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nkspin = Nkpt*Nspin

    evars = ElecVars(Ham)
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

    #println(Ham.energies)
    println("Etot = ", Etot)

    #println("rotPrevCinv:")
    #print_vec_mat(rotPrevCinv)

    #println("g Haux")
    #print_vec_mat(g.Haux)

    #d = deepcopy(Kg)
    d = deepcopy(g)

    #d.Haux[1][2,2] = 99
    #print_vec_mat(d.Haux[1:1])

    #println("Kg Haux")
    #print_vec_mat(Kg.Haux[1:1])

    # Constrain
    constrain_search_dir!( d, evars )

    β = 0.0
    gPrevUsed = true
    gKnormPrev = 0.0
    force_grad_dir = true

    # Begin iter
    NiterMax = 5
    Etot_old = Etot
    for iter in 1:NiterMax

        #gKnorm = dot_ElecGradient(g, Kg)
        #gKnorm = dot_ElecGradient(g, g)
        
        gKnorm = real(dot(g.Haux, g.Haux))
        println("gKnorm = ", gKnorm)

        gKnorm = 2*real(dot(g.psiks, Kg.psiks))
        println("gKnorm = ", gKnorm)

        #if !force_grad_dir
        #    
        #    dotgd = dot_ElecGradient(g, d)
        #    
        #    if gPrevUsed
        #        dotgPrevKg = dot_ElecGradient(gPrev, Kg)
        #        #dotgPrevKg = 2*real(dot(gPrev.psiks, Kg.psiks))
        #    else
        #        dotgPrevKg = 0.0
        #    end
        #    β = (gKnorm - dotgPrevKg)/gKnormPrev # Polak-Ribiere
        #    println("β = ", β)
        #    if β < 0.0
        #        println("Resetting β")
        #        β = 0.0
        #    end
        #end

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

        #println("rotPrevCinv")
        #print_vec_mat(rotPrevCinv[1:1])

        α = linmin_grad!( Ham, evars, g, d, kT, rotPrev, rotPrevC, rotPrevCinv )
        #α = linmin_grad!( Ham, evars.psiks, g.psiks, d.psiks )
        println("α = ", α)
        #α = 3e-5

        #do_step!( α, 0.0, evars, d, rotPrev, rotPrevC, rotPrevCinv )
        #do_step!( 0.0, α, evars, d, rotPrev, rotPrevC, rotPrevCinv )
        do_step!( α, evars, d, rotPrev, rotPrevC, rotPrevCinv )

        #println("rotPrev")
        #print_vec_mat(rotPrev[1:1])
        
        #println("rotPrevCinv")
        #print_vec_mat(rotPrevCinv[1:1])

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

        
        #println("\nevars.Haux_eigs")
        #println(evars.Haux_eigs)
        #println(Ham.electrons.Focc)

        #println("evars.Hsub")
        #print_vec_mat(evars.Hsub)
        
        #println("g Haux")
        #print_vec_mat(g.Haux)
        
        #println("Kg Haux")
        #print_vec_mat(Kg.Haux)
    end

    println("Pass here")
end

main()