using PWDFT

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

        Si  0.0   0.0   0.0
        Si  0.25  0.25  0.35
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(5.431*ANG2BOHR)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)
    write_xsf( "TEMP_Si.xsf", atoms )

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
    
    # Search for unique mapping (groups of kpoints with similar symmetry)
    umap = unique(mapping)
    println("umap = ", umap)

    kcount = zeros(Int64,num_ir)
    for ik = 1:num_ir
        kcount[ik] = count( i -> ( i == umap[ik] ), mapping )
    end

    list_ir_k = []
    for ikk = 1:num_ir
        for ik = 1:Nkpt
            if umap[ikk] == mapping[ik]
                append!( list_ir_k, [kgrid[:,ik]] )
                break
            end
        end
    end
    println( list_ir_k )

    kred = zeros(Float64,3,num_ir)
    for ik = 1:num_ir
        kred[1,ik] = list_ir_k[ik][1] / mesh[1]
        kred[2,ik] = list_ir_k[ik][2] / mesh[2]
        kred[3,ik] = list_ir_k[ik][3] / mesh[3]
    end
    println(RecVecs)
    kred = RecVecs'*kred  # convert to cartesian
    println(kred)

    # calculate the weights
    wk = kcount[:]/sum(kcount)*2.0  # spin-degenerate
    for ik = 1:num_ir
        @printf("%5d map=%d (%8.5f,%8.5f,%8.5f) %8.5f\n", ik, umap[ik],
                kred[1,ik], kred[2,ik], kred[3,ik], wk[ik])
    end

end


function test_Si_fcc_kred()
    
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        2

        Si  0.0   0.0   0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(5.431*ANG2BOHR)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)
    write_xsf( "TEMP_Si.xsf", atoms )

    meshk = [3,3,3]
    shiftk = [0,0,0]
    kred, wk = gen_kgrid_reduced( atoms, meshk, shiftk )

end

#test_Si_fcc_kred()
test_Si_fcc_get_ir()
#test_MonkhorstPack( mesh )


