using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("../create_Ham.jl")
#include("subspace_rotation.jl")
include("ElectronicVars.jl")
include("eval_L_tilde.jl")

function calc_alpha_CG!(
    g::ElectronicVars, gt::ElectronicVars, d::ElectronicVars,
    α_t::Float64,
    α_ψ::Vector{Float64},
    α_η::Vector{Float64}
)

    Nkspin = length(g.ψ)

    for i in 1:Nkspin

        denum_ψ = real(sum(conj(g.ψ[i] - gt.ψ[i]).*d.ψ[i]))
        #println("denum_ψ = ", denum_ψ)
        #println("num_ψ   = ", real(sum(conj(g.ψ[i]).*d.ψ[i])) )
        #if abs(denum) <= 1e-6
        if denum_ψ != 0.0
            #α_ψ[i] = abs( α_t*real(sum(conj(g.ψ[i]).*d.ψ[i])) / denum_ψ )
            α_ψ[i] = α_t*real(sum(conj(g.ψ[i]).*d.ψ[i])) / denum_ψ
        else
            α_ψ[i] = 0.0
        end

        denum_η = real(sum(conj(g.η[i] - gt.η[i]).*d.η[i]))
        #println("denum_η = ", denum_η)
        #println("num_η   = ", real(sum(conj(g.η[i]).*d.η[i])) )
        #if abs(denum) <= 1e-6
        if denum_η != 0.0
            #α_η[i] = abs( α_t*real(sum(conj(g.η[i]).*d.η[i])) / denum_η )
            α_η[i] = α_t*real(sum(conj(g.η[i]).*d.η[i])) / denum_η
        else
            α_η[i] = 0.0
        end

    end
    return
end

function trial_evars!( ec, e, d, α_t )
    Nkspin = length(e.ψ)
    for i in 1:Nkspin
        ec.ψ[i] = e.ψ[i] + α_t*d.ψ[i]
        ec.η[i] = e.η[i] + α_t*d.η[i]
        ec.η[i] = 0.5*( ec.η[i] + ec.η[i]' )
    end
    return
end

"""
evars_t: memory for trial evars
"""
function quadratic_line_minimizer!(Ham, evars, evars_t, d, E, g, α_t, α_1, α_2)

    Nkspin = length(evars.ψ)
    Δ = zeros(Nkspin)

    # compute directional derivative
    d_ψ, d_η = dot(d, g)
    for i in 1:Nkspin
        Δ[i] = 2.0*( d_ψ[i] + d_η[i] )
    end

    println("quadratic_line_minimizer: α_t = ", α_t)
    if sum(Δ) > 0
        println("Wrong direction !!!!")
        exit()
    end

    evars_t = copy(evars)

    trial_evars!( evars_t, evars, d, α_t )

    constraint!( Ham, evars_t )
    E_trial = eval_L_tilde!( Ham, evars_t )

    ss = sum(Δ)
    curvature = ( E_trial - ( E + α_t*ss ) ) /α_t^2
    #println("curvature = ", curvature, " num = ", num)
    cc = -ss/(2*curvature)
    if cc < 0
        #println("Wrong sign of α")
        α_1[:] .= α_t
        α_2[:] .= α_t
        #α_t = 3.0*α_t
    else
        α_1[:] .= cc
        α_2[:] .= cc
    end

    return α_t
end



function test_linmin()
    Random.seed!(1234)

    #Ham = create_Ham_atom_Pt_smearing(a=10.0)
    Ham = create_Ham_Al_fcc_smearing( meshk=[1,1,1], Nspin=1 )
    #Ham = create_Ham_Pt_fcc_smearing( meshk=[1,1,1] )

    println(Ham)

    evars = rand_ElectronicVars(Ham)

    g_evars = copy(evars)

    Ham.energies.NN = calc_E_NN( Ham.atoms )
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )

    guess_evars!( Ham, evars, NiterMax=3 )

    constraint!( Ham, evars )
    Etot_old = eval_L_tilde!( Ham, evars )

    println("\n U real part = ")
    display(real(evars.U[1]))
    println("\n U imag part = ")
    display(imag(evars.U[1]))

    #α_t = 3e-5
    α_t = 1.0

    Nkspin = length(evars.ψ)
    α_ψ = zeros(Nkspin)
    α_η = zeros(Nkspin)

    evarsc = copy(evars)
    d_evars = copy(evars)
    gt_evars = copy(evars)
    Kg_evars = copy(evars)
    d_old_evars = copy(evars)

    v1_ψ = zeros(Nkspin)
    v2_ψ = zeros(Nkspin)
    v3_ψ = zeros(Nkspin)

    v1_η = zeros(Nkspin)
    v2_η = zeros(Nkspin)
    v3_η = zeros(Nkspin)

    Nconverges = 0

    check_Hsub( Ham, evars )

    κ_0 = 1.0
    κ = κ_0

    evars_old = copy(evars)

    for iter = 1:50

        #rotate_evars!( Ham, evars )
        grad_eval_L_tilde!( Ham, evars, g_evars )

        if iter > 1
            v1_ψ, v1_η = dot( g_evars, d_evars )
            v2_ψ, v2_η = dot( g_evars, g_evars )
            v3_ψ, v3_η = dot( d_evars, d_evars )
            for i in 1:Nkspin
                @printf("cosine angle psiks = %18.10f, ", v1_ψ[i]/sqrt(v2_ψ[i]*v3_ψ[i]))
                @printf("cosine angle Haux  = %18.10f\n", v1_η[i]/sqrt(v2_η[i]*v3_η[i]))
            end
        end
        #print_Haux( evars, "evars after eval_L_tilde")
        print_Haux( g_evars, "g_evars after eval_L_tilde")

        calc_primary_search_dirs!( Ham, evars, Kg_evars, κ=κ )
        d_evars = copy( Kg_evars )

        # use unpreconditioned gradients
        # try different sign of g_evars.η
        #d_evars = ElectronicVars( -g_evars.ψ, -g_evars.η )

        #trial_evars!( evarsc, evars, d_evars, α_t )
        #constraint!( Ham, evarsc )
        #grad_eval_L_tilde!( Ham, evarsc, gt_evars )
        #print_Haux(gt_evars, "gt_evars")

        #calc_alpha_CG!( g_evars, gt_evars, d_evars, α_t, α_ψ, α_η )
        #println("iter = ", iter, " α_ψ = ", α_ψ)
        #println("iter = ", iter, " α_η = ", α_η)

        α_t = quadratic_line_minimizer!(Ham, evars, evarsc, d_evars, Etot_old, g_evars, α_t, α_ψ, α_η)

        # update evars
        axpy!( α_ψ, α_η, evars, d_evars )

        constraint!( Ham, evars )

        println("\n U real part = ")
        display(real(evars.U[1]))
        println("\n U imag part = ")
        display(imag(evars.U[1]))
        println("\n")

        Etot = eval_L_tilde!( Ham, evars )
        #print_Haux(evars, "evars after eval_L_tilde!")
        dEtot = Etot_old - Etot

        _, dd_η = dot(g_evars, g_evars)
        if sum(dd_η) > 1e-6
            if dEtot < 0.0
                evars = copy(evars_old)
                κ = 0.5*κ
                @printf("Iteration %8d Undoing step\n", iter)
                @printf("Iteration %8d Reducing κ to %18.10f\n", iter, κ)
                continue
            end
        end

        @printf("Iteration %8d %18.10f %18.10e\n", iter, Etot, Etot_old - Etot)
        println("Iteration norm grad: ", dot(g_evars, g_evars))

        check_Hsub( Ham, evars )

        if abs(Etot_old - Etot) < 1e-6
            Nconverges = Nconverges + 1
        else
            Nconverges = 0
        end
        if Nconverges >= 2
            @printf("\nEmin_Haux_linmin is converged in iter: %d\n", iter)
            break
        end

        Etot_old = Etot
        d_old_evars = copy(d_evars) # not used
        evars_old = copy(evars)
    end

    print_ebands( Ham )
    println( Ham.energies )
end
@time test_linmin()
