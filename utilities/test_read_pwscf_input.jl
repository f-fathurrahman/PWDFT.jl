using Printf
using OffsetArrays
using LinearAlgebra

using Random
Random.seed!(1234)

using PWDFT

include("PWSCFInput.jl")
pwinput = PWSCFInput(ARGS[1])
