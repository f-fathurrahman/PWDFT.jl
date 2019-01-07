using PWDFT

function main()
    atoms = Atoms(ext_xyz_file="TEMP_slab.xyz")
    println(atoms)
    write_xsf("TEMP_slab.xsf", atoms)
end

main()