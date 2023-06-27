using PWDFT
using Random
using LinearAlgebra
using Test

const DIR_PWDFT = joinpath( dirname(pathof(PWDFT)), "..")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")
const DIR_PSP_GTH_LDA = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("create_objects.jl")

@assert length(ARGS) == 1

include(ARGS[1])