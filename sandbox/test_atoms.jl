if VERSION <= v"0.6.3"
    using PWDFT
else
    push!(LOAD_PATH, pwd())
    using PWDFT07
end

function test_main()
    atoms = init_atoms_xyz("../structures/CuSO4.xyz", verbose=true)
    println(atoms)
    
    dummy_atoms = Atoms()
    println(dummy_atoms)
end

test_main()
