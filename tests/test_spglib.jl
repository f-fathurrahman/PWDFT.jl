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
   atoms = init_atoms_xyz("../structures/corrundum.xyz", in_bohr=true)
   atoms.LatVecs = zeros(3,3)
   atoms.LatVecs[:,1] = [4.8076344022756095, -2.4038172011378047, 0]
   atoms.LatVecs[:,2] = [0.0, 4.1635335244786962, 0.0]
   atoms.LatVecs[:,3] = [0.0, 0.0, 13.1172699198127543]
   atoms.LatVecs = atoms.LatVecs*ANG2BOHR
   # convert back to "true" bohr
   atoms.positions = atoms.LatVecs*atoms.positions
   println(atoms)

   write_xsf( "TEMP.xsf", atoms.LatVecs/ANG2BOHR, atoms.positions/ANG2BOHR, atsymbs=atoms.atsymbs )
   
   Natoms_prim = spg_find_primitive(atoms)
   println("Natoms_prim = ", Natoms_prim)

end

#test_BCC()
test_corrundum()


