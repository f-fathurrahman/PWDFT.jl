using Printf
using PWDFT

function test_Si_diamond()
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        2

        Si  0.0   0.0   0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(5.431*ANG2BOHR)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)

    kpoints = KPoints( atoms, [3,3,3], [0,0,0] )
    println(kpoints)
end

function test_Fe_bcc()
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        1

        Fe  0.0   0.0   0.0
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_bcc(5.0*ANG2BOHR)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)

    kpoints = KPoints( atoms, [3,3,3], [0,0,0] )
    println(kpoints)
end

function test_kpath()
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        2

        Si  0.0   0.0   0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(5.431*ANG2BOHR)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)

    kdict = get_special_kpoints("fcc")
    println(kdict)
end


function test_kpoints_from_file()
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        2

        Si  0.0   0.0   0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(5.431*ANG2BOHR)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)
    #
    kpoints = kpoints_from_file(atoms, "KPATH_FCC_60")
    println(kpoints)
    #
    pw = PWGrid(15.0, atoms.LatVecs, kpoints=kpoints)
    println(pw)
end

test_kpoints_from_file()

#test_Si_diamond()

#test_Fe_bcc()

#test_kpath()
