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

    evars = rand_ElectronicVars(Ham)

    g_evars = copy(evars)
    Kg_evars = copy(evars)

    Ham.energies.NN = calc_E_NN( Ham.atoms )
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )

    print_Haux(evars, "evars before eval_L_tilde")

    constraint!( Ham, evars )

    Etot = eval_L_tilde!( Ham, evars )

    print_Haux( evars, "evars after eval_L_tilde" )

    println("Etot = ", Etot)

    grad_eval_L_tilde!(Ham, evars, g_evars)

    print_Haux( g_evars, "g_evars before grad_eval_L_tilde" )

    calc_primary_search_dirs!( Ham, evars, Kg_evars, κ=1.0 )
    print_Haux( Kg_evars, "Kg_evars before grad_eval_L_tilde" )

    Etot = eval_L_tilde!( Ham, evars )

    println("Etot = ", Etot)
end

test_grad_eval_L_tilde()
