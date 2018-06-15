push!(LOAD_PATH, "../src")
using PWDFT

function test_main()
    atoms = init_atoms_xyz("../structures/CuSO4.xyz", verbose=true)
    println(atoms)
    
    dummy_atoms = Atoms()
    println(dummy_atoms)
end

test_main()
