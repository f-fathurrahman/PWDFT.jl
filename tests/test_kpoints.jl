using PWDFT

function test_main()
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        2

        Si  0.0   0.0   0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(5.431*ANG2BOHR)
    atoms.positions = atoms.LatVecs*atoms.positions

    kpoints = KPoints( atoms, [3,3,3], [0,0,0] )
    println(kpoints)
end

test_main()

