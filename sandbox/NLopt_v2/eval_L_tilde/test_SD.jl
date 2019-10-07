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

function test_SD()
    Random.seed!(1234)

    #Ham = create_Ham_atom_Pt_smearing(a=10.0)
    #Ham = create_Ham_Al_fcc_smearing( meshk=[1,1,1], Nspin=1, ecutwfc=15.0 )
    Ham = create_Ham_Pt_fcc_smearing( meshk=[1,1,1] )

    println(Ham)

    evars = rand_ElectronicVars(Ham)

    g_evars = copy(evars)

    Ham.energies.NN = calc_E_NN( Ham.atoms )
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )

    guess_evars!( Ham, evars, NiterMax=20 )

    constraint!( Ham, evars )
    Etot_old = eval_L_tilde!( Ham, evars )

    α_t = 3e-5
    β_t = 3e-5  # also set this?

    Nconverges = 0
    for iter = 1:10

        constraint!( Ham, evars )
        grad_eval_L_tilde!( Ham, evars, g_evars )
        print_Haux(g_evars, "g_evars after grad_eval_L_tilde")

        axpy!( -α_t, -α_t, evars, g_evars )  # XXX: try different sign for Haux?
        #axpy!( -α_t, -β_t, evars, g_evars )  # XXX: try different sign for Haux?
        #axpy!( -α_t, 1e-8, evars, g_evars )  # XXX: try different sign for Haux?

        constraint!( Ham, evars )
        Etot = eval_L_tilde!( Ham, evars )
        print_Haux(evars, "evars after grad_eval_L_tilde")

        @printf("Iteration %8d %18.10f %18.10e\n", iter, Etot, Etot_old - Etot)

        check_Hsub( Ham, evars )

        if abs(Etot_old - Etot) < 1e-10
            Nconverges = Nconverges + 1
        else
            Nconverges = 0
        end
        if Nconverges >= 3
            @printf("\nEmin_Haux_SD is converged in iter: %d\n", iter)
            break
        end

        Etot_old = Etot
    end

    print_ebands( Ham )
    println( Ham.energies )
end
@time test_SD()
