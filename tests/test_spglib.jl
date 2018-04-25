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
   
   println( atoms.LatVecs )

   write_xsf( "TEMP.xsf", atoms )
   
   Natoms_prim = spg_find_primitive( atoms )
   println( "Natoms_prim = ", Natoms_prim )

   println( atoms.LatVecs )

   red_atoms = reduce_atoms( atoms )
   println( red_atoms )
   write_xsf( "TEMP_reduced.xsf", red_atoms )

end


function test_MonkhorstPack( mesh::Array{Int64,1} )
    ik = 0
    kpts = zeros(Float64,3,prod(mesh))
    for k = 1:mesh[3]
    for j = 1:mesh[2]
    for i = 1:mesh[1]
        ik = ik + 1
        kpts[1,ik] = (2*i - mesh[1] - 1)/(2*mesh[1])
        kpts[2,ik] = (2*j - mesh[2] - 1)/(2*mesh[2])
        kpts[3,ik] = (2*k - mesh[3] - 1)/(2*mesh[3])
        @printf("%5d %8.5f %8.5f %8.5f\n", ik, kpts[1,ik], kpts[2,ik], kpts[3,ik])
    end
    end
    end
end

function test_Si_fcc_get_ir()
    
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        2

        Si  0.0  0.0  0.0
        Si  0.5  0.5  0.5
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(5.431*ANG2BOHR)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)

    RecVecs = 2*pi*inv(atoms.LatVecs')/(2*pi/(5.431*ANG2BOHR))
    b1 = RecVecs[:,1]
    b2 = RecVecs[:,2]
    b3 = RecVecs[:,3]

    mesh = [3,3,3]
    is_shift = [0,0,0]
    num_ir, kgrid, mapping =
    spg_get_ir_reciprocal_mesh( atoms, mesh, is_shift, is_time_reversal=1 )

    Nkpt = prod(mesh)
    println("num_ir = ", num_ir)
    for ik = 1:Nkpt
        kpt = b1*kgrid[1,ik]/mesh[1] + b2*kgrid[2,ik]/mesh[2] + b3*kgrid[3,ik]/mesh[3]
        @printf("%5d: %8.5f %8.5f %8.5f  => %d\n", ik, kpt[1], kpt[2], kpt[3], mapping[ik])
#        @printf("%5d: %8.5f %8.5f %8.5f  => %d\n", ik,
#                kgrid[1,ik]/mesh[1], kgrid[2,ik]/mesh[2], kgrid[3,ik]/mesh[3], mapping[ik])
    end
    println("kgrid = ", size(kgrid))
    println("mapping = ", size(mapping))

    umap = unique(mapping)
    kcount = zeros(Int64,num_ir)
    for ik = 1:num_ir
        kcount[ik] = count( i -> ( i == umap[ik] ), mapping )
    end

    println(umap)
    # calculate the weights
    wk = kcount[:]/sum(kcount)*2.0  # spin-degenerate
    for ik = 1:num_ir
        @printf("%5d map=%d %8.5f\n", ik, umap[ik], wk[ik])
    end

    test_MonkhorstPack( mesh )

end

test_Si_fcc_get_ir()
#test_BCC()
#test_corrundum()


