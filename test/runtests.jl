using PWDFT

using LinearAlgebra
using Test

include("test_gen_lattice.jl")
include("test_atoms.jl")

@test gen_lattice_cubic(16.0) == 16.0*Matrix(Diagonal(ones(3)))

@test test_gen_lattice() == nothing

@test test_atoms() == nothing

println("")
println("All test successfully run")
println("")

