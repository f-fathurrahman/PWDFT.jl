using PWDFT

function test01()
    atoms = init_atoms_xyz("../structures/CuSO4.xyz", verbose=true)
    println(atoms)
    
    dummy_atoms = Atoms()
    println(dummy_atoms)
end

function test02()
    atoms = Atoms(xyz_file="../structures/CuSO4.xyz")
    println(atoms)
end


test02()
