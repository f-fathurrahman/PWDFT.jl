using Printf
using LinearAlgebra

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("SymmetryBase.jl")

function create_Si_fcc()
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(5.431*ANG2BOHR))
    # FIXME: Zvals is not set
    return atoms
end

function create_Si_fcc_cubic()
    atoms = Atoms(xyz_string_frac=
        """
        8

        Si   0.000000000000000   0.000000000000000   0.000000000000000 
        Si   0.750000000000000   0.750000000000000   0.250000000000000 
        Si   0.500000000000000   0.000000000000000   0.500000000000000 
        Si   0.750000000000000   0.250000000000000   0.750000000000000 
        Si   0.000000000000000   0.500000000000000   0.500000000000000 
        Si   0.250000000000000   0.250000000000000   0.250000000000000 
        Si   0.250000000000000   0.750000000000000   0.750000000000000 
        Si   0.500000000000000   0.500000000000000   0.000000000000000 
        """, in_bohr=true, LatVecs=gen_lattice_cubic(5.43070*ANG2BOHR))
    return atoms
end

function main()
    
    #atoms = create_Si_fcc()
    atoms = create_Si_fcc_cubic()
    
    sym_base = SymmetryBase(atoms)
    @time sym_base = SymmetryBase(atoms)

    #println(sym_base)

    Natoms = atoms.Natoms
    v = zeros(Float64,3,Natoms)
    v[:,1] = [1.0, 0.0, 0.0]
    v[:,2] = [0.0, 1.0, 1.0]

    symmetrize_vector!( atoms.LatVecs, sym_base, v )
    println("v = ")
    display(v); println()
end

main()