using Printf
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
    println(atoms)

    ecutwfc_Ry = 30.0
    kpoints = KPoints( atoms, [4,4,4], [0,0,0] )
    pw = PWGrid( ecutwfc_Ry*0.5, atoms.LatVecs, kpoints=kpoints )
    println(pw)
end

@time test_main()
