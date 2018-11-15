import InteractiveUtils
using Printf
using Random
using BenchmarkTools
using PWDFT

include("time_gen_lattice.jl")
include("time_Atoms.jl")
include("time_PWGrid.jl")
include("time_Hamiltonian.jl")
include("time_op_kpt_1.jl")

println()
InteractiveUtils.versioninfo()

#time_gen_lattice()
#time_Atoms()
#time_PWGrid()
#time_Hamiltonian()
time_op_kpt_1()
