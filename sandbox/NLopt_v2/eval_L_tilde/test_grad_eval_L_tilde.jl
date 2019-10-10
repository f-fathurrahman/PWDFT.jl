using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("../create_Ham.jl")
include("ElectronicVars.jl")
include("eval_L_tilde.jl")

function test_grad_eval_L_tilde()

    Random.seed!(1234)

    Ham = create_Ham_Al_fcc_smearing( meshk=[1,1,1] )
    #Ham = create_Ham_Pt_fcc_smearing( meshk=[1,1,1] )

    evars = rand_ElectronicVars(Ham)

    g_evars = copy(evars)
    Kg_evars = copy(evars)

    Ham.energies.NN = calc_E_NN( Ham.atoms )
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )

    #guess_evars!( Ham, evars, NiterMax=20 )

    print_Haux(evars, "evars before eval_L_tilde")

    constraint!( Ham, evars )

    println("\n U real part = ")
    display(real(evars.U[1]))
    println("\n U imag part = ")
    display(imag(evars.U[1]))

    Etot = eval_L_tilde!( Ham, evars )
    print_Haux( evars, "evars after eval_L_tilde" )

    println("Etot = ", Etot)

    grad_eval_L_tilde!(Ham, evars, g_evars)
    print_Haux( g_evars, "g_evars after grad_eval_L_tilde" )

    calc_primary_search_dirs!( Ham, evars, Kg_evars, κ=1.0 )
    print_Haux( Kg_evars, "Kg_evars after calc_primary_search_dirs" )

    Etot = eval_L_tilde!( Ham, evars )
    println("Etot = ", Etot)

    check_Hsub( Ham, evars )

    print_ebands(Ham)

    dd_ψ, dd_η = dot( evars, g_evars )
    println("dot_ψ = ", dd_ψ)
    println("dot_η = ", dd_η)

    dd_ψ, dd_η = dot( g_evars, g_evars )
    println("norm_g_ψ = ", dd_ψ)
    println("norm_g_η = ", dd_η)

    dd_ψ, dd_η = dot( evars, Kg_evars )
    println("prec dot_ψ = ", dd_ψ)
    println("prec dot_η = ", dd_η)

end

test_grad_eval_L_tilde()
