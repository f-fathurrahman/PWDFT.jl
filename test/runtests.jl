using PWDFT

using LinearAlgebra
using Test

include("test_gen_lattice.jl")
include("test_atoms.jl")

@test gen_lattice_cubic(16.0) == 16.0*Matrix(Diagonal(ones(3)))

@test test_gen_lattice() == nothing

@test test_atoms() == nothing

include("test_Hamiltonian.jl")

@test test_Hamiltonian_v1() == nothing
@test test_Hamiltonian_v2() == nothing

include("test_SCF.jl")

@test test_SCF_v1() == nothing
@test test_SCF_v2() == nothing

println("")
println("All test successfully run")
println("")

