using Printf
using LinearAlgebra
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

function test_Si_fcc()
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0   0.0   0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true,
        LatVecs=gen_lattice_fcc(5.431*ANG2BOHR))
    
    println(atoms)

    println(SymmetryInfo(atoms))
end


function test_GaAs()
    atoms = Atoms(xyz_string_frac=
        """
        2

        Ga  0.0   0.0   0.0
        As  0.25  0.25  0.25
        """, in_bohr=true,
        LatVecs=gen_lattice_fcc(10.6839444516))
    
    println(atoms)
    println(SymmetryInfo(atoms))
end

function test_H_fcc()
    atoms = Atoms(xyz_string_frac=
        """
        1

        H  0.0   0.0   0.0
        """, in_bohr=true,
        LatVecs=gen_lattice_fcc(5.0))
    
    println(atoms)
    println(SymmetryInfo(atoms))
end


function test_CH4()
    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "CH4.xyz"),
                   LatVecs=gen_lattice_cubic(16.0))
    println(atoms)
    println(SymmetryInfo(atoms))
end


test_Si_fcc()
#test_GaAs()
#test_CH4()
#test_H_fcc()