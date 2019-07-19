using PWDFT
using Random
using LinearAlgebra
using Test

const DIR_PWDFT = joinpath( dirname(pathof(PWDFT)), "..")

include("test_gen_lattice.jl")
include("test_atoms.jl")
include("test_Hamiltonian.jl")

include("test_xc.jl")
include("test_spglib.jl")

include("H_atom.jl")
include("Si_fcc.jl")
include("Fe_bcc_PBE_spinpol.jl")

println("")
println("All test successfully run")
println("")
