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
    println(atoms)

    Natoms_prim = spg_find_primitive(atoms)
    println("Natoms_prim = ", Natoms_prim)

    println(atoms)
end


function test_corrundum()

   atoms = init_atoms_xyz("../structures/corrundum.xyz")
   
   atoms.LatVecs = zeros(3,3)
   atoms.LatVecs[:,1] = [4.1228005, -2.3803000, 0.0000000]
   atoms.LatVecs[:,2] = [0.0000000, 4.7606000, 0.0000000]
   atoms.LatVecs[:,3] = [0.0000000, 0.0000000, 12.9940000]
   atoms.LatVecs = atoms.LatVecs*ANG2BOHR
   
   println(atoms.LatVecs)

   write_xsf( "TEMP.xsf", atoms )
   
   Natoms_prim = spg_find_primitive(atoms)
   println("Natoms_prim = ", Natoms_prim)

   println(atoms.LatVecs)

   red_atoms = reduce_atoms(atoms)
   println(red_atoms)
   write_xsf( "TEMP_reduced.xsf", red_atoms )

end

#test_BCC()
test_corrundum()


