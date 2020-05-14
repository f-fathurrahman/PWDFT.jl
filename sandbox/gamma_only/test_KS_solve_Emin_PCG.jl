using Printf
using LinearAlgebra
using Random
using FFTW

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("PWGridGamma.jl")
include("wrappers_fft_gamma.jl")
include("ortho_GS_gamma.jl")
include("PsPotNLGamma.jl")
include("HamiltonianGamma.jl")
include("BlochWavefuncGamma.jl")
include("calc_rhoe_gamma.jl")
include("Poisson_solve_gamma.jl")
include("op_K_gamma.jl")
include("op_V_loc_gamma.jl")
include("op_V_Ps_nloc_gamma.jl")
include("op_H_gamma.jl")
include("calc_energies_gamma.jl")
include("calc_grad_gamma.jl")

include("setup_guess_wavefunc.jl")

include("KS_solve_Emin_PCG_dot.jl")
include("calc_energies_grad.jl")
include("linmin_grad.jl")

include("KS_solve_Emin_PCG_dot_gamma.jl")
include("calc_energies_grad_gamma.jl")
include("linmin_grad_gamma.jl")

include("unfold_BlochWavefuncGamma.jl")

function test_01()

    Random.seed!(1234)

    #atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "H2.xyz"),
    #               LatVecs = gen_lattice_sc(16.0) )
    #pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    
    atoms = Atoms( ext_xyz_file=joinpath(DIR_STRUCTURES, "NH3.xyz") )
    pspfiles = [joinpath(DIR_PSP, "N-q5.gth"),
                joinpath(DIR_PSP, "H-q1.gth")]

    # Initialize Hamiltonian
    ecutwfc = 15.0
    Ham = HamiltonianGamma( atoms, pspfiles, ecutwfc )

    Ham_ = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )

    psis = randn_BlochWavefuncGamma(Ham)
    psiks = unfold_BlochWavefuncGamma( Ham.pw, Ham_.pw, psis )

    #KS_solve_Emin_PCG_dot!( Ham, psis, NiterMax=50 )

    KS_solve_Emin_PCG_dot!( Ham_, psiks, startingrhoe=:random, skip_initial_diag=true, NiterMax=50 )
    #KS_solve_Emin_PCG!( Ham_, psiks, startingrhoe=:random, skip_initial_diag=true )

end

#=
Kinetic    energy:       8.0785591277
Ps_loc     energy:     -31.8301682150
Ps_nloc    energy:       0.9804183905
Hartree    energy:      11.8324308638
XC         energy:      -3.4421445627
-------------------------------------
Electronic energy:     -14.3809043957
NN         energy:       3.1331525188
-------------------------------------
Total      energy:     -11.2477518770
=#

test_01()