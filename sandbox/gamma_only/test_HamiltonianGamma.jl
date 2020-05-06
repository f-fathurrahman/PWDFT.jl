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
include("BlochWavefuncGamma.jl")
include("PsPotNLGamma.jl")
include("HamiltonianGamma.jl")

function test_01()
    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "H2.xyz"),
                   LatVecs = gen_lattice_sc(16.0) )
    
    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0
    Ham = HamiltonianGamma( atoms, pspfiles, ecutwfc )
    
    Ham_ = Hamiltonian( atoms, pspfiles, ecutwfc )

    println(sum(Ham.potentials.Ps_loc))

    println(sum(Ham_.potentials.Ps_loc))

    println("Pass here")
end

test_01()