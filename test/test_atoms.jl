function test_atoms()
    dir_PWDFT = joinpath( dirname(pathof(PWDFT)), "..")
    atoms = init_atoms_xyz(
        joinpath( dir_PWDFT, "structures/CuSO4.xyz")
    )
    atoms.LatVecs = gen_lattice_cubic(10.0)
    println(atoms)

    atoms = init_atoms_xyz_string(
        """
        1

        H  0.0  0.0  0.0
        """
    )
    atoms.LatVecs = gen_lattice_cubic(10.0)
    println(atoms)

    return nothing
end