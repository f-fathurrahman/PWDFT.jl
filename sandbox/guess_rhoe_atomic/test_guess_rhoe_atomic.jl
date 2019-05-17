using Printf
using Random
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("guess_rhoe_atomic.jl")

function init_Ham_GaAs()
    atoms = Atoms(xyz_string_frac=
        """
        2

        Ga  0.0  0.0  0.0
        As  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.6839444516))

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Ga-q3.gth"),
                joinpath(DIR_PSP, "As-q5.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
end

function main()
    Ham = init_Ham_GaAs()
    Rhoe1 = guess_rhoe_atomic( Ham )
    Rhoe2 = guess_rhoe( Ham )

    println("diff Rhoe = ", sum(abs.(Rhoe1[:,1] - Rhoe2)))
end

main()
