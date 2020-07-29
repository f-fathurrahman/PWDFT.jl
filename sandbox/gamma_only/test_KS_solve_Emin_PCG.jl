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

    #atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "CO2.xyz"), LatVecs=gen_lattice_sc(16.0) )
    #pspfiles = [joinpath(DIR_PSP, "C-q4.gth"),
    #            joinpath(DIR_PSP, "O-q6.gth")]

    #atoms = Atoms( ext_xyz_file=joinpath(DIR_STRUCTURES, "H2O.xyz") )
    #pspfiles = [joinpath(DIR_PSP, "O-q6.gth"),
    #            joinpath(DIR_PSP, "H-q1.gth")]

    # Initialize Hamiltonian
    ecutwfc = 15.0
    Ham = HamiltonianGamma( atoms, pspfiles, ecutwfc )

    Ham_ = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )

    psis = randn_BlochWavefuncGamma(Ham)
    psiks = unfold_BlochWavefuncGamma( Ham.pw, Ham_.pw, psis )

    @time KS_solve_Emin_PCG_dot!( Ham, psis, NiterMax=200 )

    #@time KS_solve_Emin_PCG_dot!( Ham_, psiks, startingrhoe=:random, skip_initial_diag=true, NiterMax=50 )
    #@time KS_solve_Emin_PCG!( Ham_, psiks, startingrhoe=:random, skip_initial_diag=true, restrict_linmin=true )
    #@time KS_solve_Emin_PCG!( Ham_, psiks, restrict_linmin=true )

end

#=

NH3:

-------------------------
Final Kohn-Sham energies:
-------------------------

Kinetic    energy:       8.0684717780
Ps_loc     energy:     -31.8615867485
Ps_nloc    energy:       0.9939824683
Hartree    energy:      11.8381690079
XC         energy:      -3.4429618488
-------------------------------------
Electronic energy:     -14.4039253430
NN         energy:       3.1331525188
-------------------------------------
Total      energy:     -11.2707728243
=#

test_01()