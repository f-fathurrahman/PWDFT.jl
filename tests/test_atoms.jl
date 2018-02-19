using PWDFT

function test_main()
    atoms = init_atoms_xyz("CuSO4.xyz", verbose=true)
    println(atoms)
end

test_main()
