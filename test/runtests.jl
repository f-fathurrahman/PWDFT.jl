using PWDFT

using LinearAlgebra
using Test

const DIR_PWDFT = joinpath( dirname(pathof(PWDFT)), "..")


#include("test_gen_lattice.jl")
#include("test_atoms.jl")
#include("test_Hamiltonian.jl")

#include("test_xc.jl")
#include("test_spglib.jl")

#include("test_KS_solve_SCF.jl")

#include("test_KS_solve_SCF_spin.jl")

#include("H_atom.jl")
include("Si_fcc.jl")

println("")
println("All test successfully run")
println("")
