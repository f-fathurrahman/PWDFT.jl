using Random
using Printf
using PWDFT

Random.seed!(1234)

include("my_scf.jl")

const DIR_PWDFT = joinpath( dirname(pathof(PWDFT)), "..")

function test_my_scf(β::Float64)
    atoms = Atoms(
        xyz_file=joinpath(DIR_PWDFT, "structures/H2.xyz"),
        LatVecs=gen_lattice_sc(16.0))
    pspfiles = [joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth", "H-q1.gth")]
    Ham = Hamiltonian( atoms, pspfiles, 15.0 )
    println(Ham)

    my_scf!(Ham, betamix=β)

    println()
    println("Kohn-Sham energy components\n")
    println(Ham.energies)
end

for β in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    test_my_scf(β)
end