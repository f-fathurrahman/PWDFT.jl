using LinearAlgebra
using Printf
using Random

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include(joinpath(DIR_PWDFT, "sandbox", "KS_solve_SCF_NLsolve.jl"))

include("MinimizeParams.jl")
include("calc_energies_grad.jl")
include("create_Ham.jl")
include("KS_solve_Emin_PCG_new.jl")
include("linmin_quad.jl")
include("linmin_grad.jl")
include("linmin_debug.jl")
include("linmin_armijo.jl")
include("KS_solve_Emin_PCG_dot.jl")

function main()
    Random.seed!(1234)

    #Ham = create_Ham_H2()
    #Ham = create_Ham_H_atom()
    #Ham = create_Ham_Si_fcc()
    Ham = create_Ham_GaAs()
    #Ham = create_Ham_NH3()
    #Ham = create_Ham_ZnO()

    psiks = rand_BlochWavefunc( Ham )

    Rhoe = calc_rhoe(Ham, psiks)
    update!(Ham, Rhoe)

    g = zeros_BlochWavefunc( Ham )
    gPrev = zeros_BlochWavefunc( Ham )
    Kg = zeros_BlochWavefunc( Ham )
    d = zeros_BlochWavefunc( Ham )
    d_old = zeros_BlochWavefunc( Ham )

    Ham.energies.NN = calc_E_NN( Ham.atoms )

    Etot = calc_energies_grad!( Ham, psiks, g, Kg )
    println("Etot = ", Etot)

    Etot_old = Etot
    Nconverges = 0
    etot_conv_thr = 1e-6
    gKnormPrev = 0.0
    β = 0.0

    Nkspin = length(psiks)
    Nstates = size(psiks[1],2)

    minim_params = MinimizeParams()
    αt = minim_params.αt_start
    α = 0.0

    for iter in 1:50

        println("\nBegin iter: ", iter)

        gKnorm = 2*real( dot(g, Kg) )
        
        if iter > 1
            #dotgd = real( dot(g, d) )
            dotgPrevKg = 2*real( dot(gPrev, Kg) )
            #println("dot g d    = ", dotgd)
            println("gKnorm     = ", gKnorm)
            println("dotgPrevKg = ", dotgPrevKg)
            β = (gKnorm - dotgPrevKg)/gKnormPrev # Polak-Ribiere
            println("β = ", β)
            if β < 0.0
                println("Resetting β")
                β = 0.0
            end
        end

        d_old = deepcopy(d)
        d = -Kg + β*d_old

        constrain_search_dir!(d, psiks)

        #_, α = linmin_grad!( Ham, psiks, g, d, Etot )
        
        linmin_success, α, αt = linmin_quad!( Ham, psiks, g, d, α, αt, Etot, minim_params )
        println("linmin_success = ", linmin_success)
        @printf("α = %e     αt = %e\n", α, αt)
        # Update
        if linmin_success
            αt = α
        end

        do_step!( psiks, α, d )

        gPrev = deepcopy(g)
        gKnormPrev = gKnorm

        Etot = calc_energies_grad!( Ham, psiks, g, Kg )

        diffE = Etot_old - Etot
        norm_grad = norm(g)/(Nkspin*Nstates)
        @printf("iter = %3d  %18.10f  %18.10e  normg = %e\n", iter, Etot, diffE, norm_grad)
        
        if abs(diffE) < etot_conv_thr
            Nconverges = Nconverges + 1
        else
            Nconverges = 0
        end

        #if (Nconverges >= 2) && (norm_grad >= 1e-5)
        #    println("Probably early convergence, continuing ...")
        #    Nconverges = 0
        #end
        
        if Nconverges >= 2
            @printf("\nConverged in iter: %d\n", iter)
            break
        end

        Etot_old = Etot
    end

end

main()