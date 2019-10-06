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
    Ham = create_Ham_Al_fcc_smearing( meshk=[1,1,1], Nspin=1, ecutwfc=30.0 )
    #Ham = create_Ham_Pt_fcc_smearing()

    println(Ham)

    evars = rand_ElectronicVars(Ham)

    g_evars = copy(evars)

    Ham.energies.NN = calc_E_NN( Ham.atoms )
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )

    # prepare guess wavefunc
    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Rhoe = zeros(Float64,Npoints,Nspin)
    if Nspin == 1
        Rhoe[:,1] = guess_rhoe( Ham )
    else
        Rhoe = guess_rhoe_atomic( Ham, starting_magnetization=[0.1] )
    end
    update!(Ham, Rhoe)

    Ham.electrons.ebands[:,:] = diag_LOBPCG!( Ham, evars.ψ, verbose=true, NiterMax=20 )
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        i = ik + (ispin -1)*Nkpt
        evars.η[i] = diagm( 0 => Ham.electrons.ebands[:,i] )
    end

    constraint!( Ham, evars )
    Etot_old = eval_L_tilde!( Ham, evars )

    α_t = 3e-5
    #β_t = 1e-1  # not good

    Nconverges = 0
    for iter = 1:10

        constraint!( Ham, evars )
        grad_eval_L_tilde!( Ham, evars, g_evars )
        print_Haux(g_evars, "g_evars after grad_eval_L_tilde")

        axpy!( -α_t, -α_t, evars, g_evars )

        constraint!( Ham, evars )
        Etot = eval_L_tilde!( Ham, evars )

        @printf("Iteration %8d %18.10f %18.10e\n", iter, Etot, Etot_old - Etot)
        if abs(Etot_old - Etot) < 1e-6
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
