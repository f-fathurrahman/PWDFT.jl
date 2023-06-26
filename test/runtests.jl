using PWDFT
using Random
using LinearAlgebra
using Test

const DIR_PWDFT = joinpath( dirname(pathof(PWDFT)), "..")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")
const DIR_PSP_GTH_LDA = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("test_gen_lattice.jl")

include("test_atoms.jl")

include("test_xc.jl")
include("test_xc_internal.jl")

include("test_Electrons.jl")

include("test_Hamiltonian.jl")

# Disabled, probaly use Spglib.jl instead of LibSymspg
#include("test_spglib.jl")

println("")
println("All test successfully run")
println("")
