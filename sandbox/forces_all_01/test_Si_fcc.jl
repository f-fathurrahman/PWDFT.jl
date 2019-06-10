using LinearAlgebra
using Printf
using Random
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT,"pseudopotentials","pade_gth")

include("calc_forces_finite_diff.jl")

function test_main()

    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))
    # LatVecs=gen_lattice_fcc(5.431*ANG2BOHR)

    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 15.0

    forces = calc_forces_finite_diff(atoms, pspfiles, ecutwfc, [3,3,3])*2.0

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