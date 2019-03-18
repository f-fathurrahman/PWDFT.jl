using LinearAlgebra
using Printf
using Random
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT,"pseudopotentials","pade_gth")

include("calc_forces_finite_diff.jl")

function test_main()

    # Atoms
    atoms = Atoms(xyz_string=
        """
        2

        H      3.83653478       4.23341768       4.23341768
        H      4.63030059       4.23341768       4.23341768
        """, LatVecs = gen_lattice_sc(16.0))

    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0

    forces = calc_forces_finite_diff(atoms, pspfiles, ecutwfc, [1,1,1])*2.0

    println("")
    println("forces (in Ry/au) = ")
    Natoms = atoms.Natoms
    atsymbs = atoms.atsymbs
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                forces[1,ia], forces[2,ia], forces[3,ia] )
    end

end

test_main()
