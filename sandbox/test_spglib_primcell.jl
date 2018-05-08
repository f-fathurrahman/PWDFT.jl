using PWDFT

function test_BCC()
    atoms = init_atoms_xyz_string(
    """
    2

    Fe  0.0   0.0   0.0
    Fe  3.0   3.0   3.0
    """
    , in_bohr=true)
    atoms.LatVecs = 6.0*eye(3)

    write_xsf( "TEMP_bcc.xsf", atoms )

    Natoms_prim = spg_find_primitive(atoms)
    println("Natoms_prim = ", Natoms_prim)

    red_atoms = reduce_atoms( atoms )
    println( red_atoms )
    write_xsf( "TEMP_bcc_reduced.xsf", red_atoms )

end


function test_Si()
    # initialize with coordinates in angstrom
    atoms = init_atoms_xyz_string("""
    8

Si   0.000000000000000   0.000000000000000   0.000000000000000 
Si   0.000000000000000   2.715350000000000   2.715350000000000 
Si   1.357675000000000   1.357675000000000   1.357675000000000 
Si   2.715350000000000   0.000000000000000   2.715350000000000 
Si   1.357675000000000   4.073024999999999   4.073024999999999 
Si   4.073024999999999   4.073024999999999   1.357675000000000 
Si   4.073024999999999   1.357675000000000   4.073024999999999 
Si   2.715350000000000   2.715350000000000   0.000000000000000
""")
   atoms.LatVecs = 5.4307*ANG2BOHR*eye(3)
   println(atoms)
   write_xsf( "TEMP_Si.xsf", atoms )

   Natoms_prim = spg_find_primitive( atoms )
   println( "Natoms_prim = ", Natoms_prim )

   red_atoms = reduce_atoms( atoms )
   println( red_atoms )
   write_xsf( "TEMP_Si_reduced.xsf", red_atoms )
end


function test_corundum()

   atoms = init_atoms_xyz("../structures/corundum.xyz")
   
   atoms.LatVecs = zeros(3,3)
   atoms.LatVecs[:,1] = [4.1228005, -2.3803000, 0.0000000]
   atoms.LatVecs[:,2] = [0.0000000, 4.7606000, 0.0000000]
   atoms.LatVecs[:,3] = [0.0000000, 0.0000000, 12.9940000]
   atoms.LatVecs = atoms.LatVecs*ANG2BOHR

   write_xsf( "TEMP_corundum.xsf", atoms )
   
   Natoms_prim = spg_find_primitive( atoms )
   println( "Natoms_prim = ", Natoms_prim )

   red_atoms = reduce_atoms( atoms )
   println( red_atoms )
   write_xsf( "TEMP_corundum_reduced.xsf", red_atoms )

end


test_Si()

#test_BCC()

#test_corundum()


