using Random
using Printf
using PWDFT

Random.seed!(1234)

include("my_scf.jl")

const DIR_PWDFT = joinpath( dirname(pathof(PWDFT)), "..")

function test_my_scf()

    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth", "Si-q4.gth")]
    ecutwfc = 20.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[4,4,4] )
    println(Ham)

    my_scf!(Ham, betamix=0.5)

    println()
    println("Kohn-Sham energy components\n")
    println(Ham.energies)
end
@time test_my_scf()
@time test_my_scf()