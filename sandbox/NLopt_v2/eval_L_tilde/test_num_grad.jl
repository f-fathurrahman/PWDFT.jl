using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("../create_Ham.jl")
include("ElectronicVars.jl")
include("eval_L_tilde.jl")

function test_num_grad()

    Random.seed!(1234)

    Ham = create_Ham_Al_fcc_smearing( meshk=[1,1,1] )
    #Ham = create_Ham_Pt_fcc_smearing( meshk=[1,1,1] )

    evars = rand_ElectronicVars(Ham)

    g_evars = copy(evars)
    Kg_evars = copy(evars)

    Ham.energies.NN = calc_E_NN( Ham.atoms )
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )

    #guess_evars!( Ham, evars, NiterMax=20 )
    evars_orig = copy(evars)

    print_Haux(evars, "evars before eval_L_tilde")

    constraint!( Ham, evars )
    #evars_orig = copy(evars) # use the diagonalized evars

    Etot = eval_L_tilde!( Ham, evars )

    #print_Haux( evars, "evars after eval_L_tilde" )
    println("Etot = ", Etot)

    Δ = 1e-5

    m = 1
    n = 2
    # vary η
    evarsc = copy(evars_orig)
    evarsc.η[1][m,n] = evars_orig.η[1][m,n] + Δ
    if m != n
        evarsc.η[1][n,m] = evars_orig.η[1][n,m] + Δ
    end
    #evarsc.η[1] = 0.5*(evarsc.η[1]' + evarsc.η[1]) # symmetrize
    print_Haux( evarsc, "evarsc plus before constraint", imag_part=false )
    constraint!(Ham, evarsc)
    print_Haux( evarsc, "evarsc plus after constraint", imag_part=false )
    Etot_plus = eval_L_tilde!( Ham, evarsc )
    println("Etot plus  = ", Etot_plus)
    print_ebands( Ham )

    evarsc = copy(evars_orig)
    evarsc.η[1][m,n] = evars_orig.η[1][m,n] - Δ
    if m != n
        evarsc.η[1][n,m] = evars_orig.η[1][n,m] - Δ
    end
    #evarsc.η[1] = 0.5*(evarsc.η[1]' + evarsc.η[1]) # symmetrize
    print_Haux( evarsc, "evarsc minus before constraint", imag_part=false )
    constraint!(Ham, evarsc)
    print_Haux( evarsc, "evarsc minus after constraint", imag_part=false )
    Etot_minus = eval_L_tilde!( Ham, evarsc )
    println("Etot minus = ", Etot_minus)
    print_ebands( Ham )

    println("deriv = ", 0.5*(Etot_plus - Etot_minus)/Δ)

    grad_eval_L_tilde!(Ham, evars, g_evars)
    print_Haux( g_evars, "g_evars after grad_eval_L_tilde" )

    #calc_primary_search_dirs!( Ham, evars, Kg_evars, κ=1.0 )
    #print_Haux( Kg_evars, "Kg_evars after calc_primary_search_dirs" )

    #Etot = eval_L_tilde!( Ham, evars )
    #println("Etot = ", Etot)

    #check_Hsub( Ham, evars )

end

test_num_grad()
