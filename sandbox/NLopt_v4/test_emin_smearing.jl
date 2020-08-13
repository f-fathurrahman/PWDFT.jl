using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("smearing.jl")
include("create_Ham.jl")
include("ElecVars.jl")
include("MinimizeParams.jl")
include("test_ElecVars.jl")
include("emin_smearing.jl")
include("linmin_grad.jl")
include("linmin_quad.jl")
include("setup_guess_wavefunc.jl")
include("KS_solve_Emin_SD_Haux.jl")
include("KS_solve_Emin_PCG_Haux.jl")
include("KS_solve_Emin_PCG_Haux_v1.jl")
include("KS_solve_Emin_PCG_Haux_v2.jl")

function main()

    Random.seed!(12345)

    kT = 0.01
    #Ham = create_Ham_atom_Si_smearing()
    #Ham = create_Ham_atom_Al_smearing()
    Ham = create_Ham_atom_C_smearing()
    #Ham = create_Ham_Al_fcc_smearing()
    #Ham = create_Ham_atom_Pt_smearing()
    #Ham = create_Ham_Pt_fcc_smearing()
    #println(Ham)

    #test_ElecVars(Ham)

    psiks = rand_BlochWavefunc(Ham)
    #setup_guess_wavefunc!( Ham, psiks, startingrhoe=:gaussian, skip_initial_diag=false )
    setup_guess_wavefunc!( Ham, psiks, startingrhoe=:random, skip_initial_diag=false )
    evars = ElecVars(Ham, psiks)
    
    #evars = ElecVars(Ham)
    
    Ham.electrons.ebands[:] = evars.Hsub_eigs[:] # Initialize Haux to eigvals of Hsub
    
    println("Initial guess:")
    print_ebands_Hsub_eigs(Ham, evars)

    KS_solve_Emin_SD_Haux!( Ham, evars, NiterMax=50 )
    #KS_solve_Emin_PCG_Haux!( Ham, evars, NiterMax=10 )
    #KS_solve_Emin_PCG_Haux_v1!( Ham, evars, NiterMax=30 )
    #KS_solve_Emin_PCG_Haux_v2!( Ham, evars, NiterMax=50 )
    #mTS = calc_entropy( wk, kT, Ham.electrons.ebands, E_fermi, Nspin )

    print_ebands(Ham.electrons, Ham.pw.gvecw.kpoints)

end

main()