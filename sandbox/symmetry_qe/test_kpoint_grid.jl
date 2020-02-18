using Printf
using LinearAlgebra

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("SymmetryBase.jl")
include("gen_kpoint_grid.jl")

function create_Si_fcc()
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.1  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(5.431*ANG2BOHR))
    # FIXME: Zvals is not set
    return atoms
end

function main()
    atoms = create_Si_fcc()
    sym_base = SymmetryBase(atoms)

    time_reversal = true

    gen_kpoint_grid( atoms.LatVecs, 3, 3, 3, 0, 0, 0, time_reversal, sym_base.s )

    kpoints = KPoints( atoms, [3,3,3], [0,0,0] )
    println(kpoints)

end

main()