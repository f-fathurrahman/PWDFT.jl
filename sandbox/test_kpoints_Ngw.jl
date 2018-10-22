using Printf
using LinearAlgebra
using PWDFT

function kpoints_Ngw( atoms:: Atoms, meshk::Array{Int64,1} )

    LatVecs = atoms.LatVecs

    RecVecs = 2*pi*inv(LatVecs')

    b1 = RecVecs[:,1]
    b2 = RecVecs[:,2]
    b3 = RecVecs[:,3]

    is_shift = [0,0,0]
    num_ir, kgrid, mapping =
    spg_get_ir_reciprocal_mesh( atoms, meshk, is_shift, is_time_reversal=1 )

    Nkpt = prod(meshk)
    for ik = 1:Nkpt
        kpt = b1*kgrid[1,ik]/meshk[1] + b2*kgrid[2,ik]/meshk[2] + b3*kgrid[3,ik]/meshk[3]
        @printf("%5d: %8.5f %8.5f %8.5f  => %d\n", ik, kpt[1], kpt[2], kpt[3], mapping[ik])
    end
    
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
                break  # from ik loop
            end
        end
    end
    println( list_ir_k )

    kred = zeros(Float64,3,num_ir)
    for ik = 1:num_ir
        kred[1,ik] = list_ir_k[ik][1] / meshk[1]
        kred[2,ik] = list_ir_k[ik][2] / meshk[2]
        kred[3,ik] = list_ir_k[ik][3] / meshk[3]
    end
    println(RecVecs)
    kred = RecVecs*kred  # convert to cartesian
    println(kred)

    # calculate the weights
    wk = kcount[:]/sum(kcount)*2.0  # spin-degenerate
    for ik = 1:num_ir
        @printf("%5d map=%d (%8.5f,%8.5f,%8.5f) %8.5f\n", ik, umap[ik],
                kred[1,ik], kred[2,ik], kred[3,ik], wk[ik])
    end

    # the group of all
    vector_all = []
    for ikr = 1:num_ir
        list_ikr = []
        for ik = 1:Nkpt
            if umap[ikr] == mapping[ik]
                append!( list_ikr, [kgrid[:,ik]] )
            end
        end
        append!(vector_all, [list_ikr])
    end
    
    RecVecs = 2*pi*inv(LatVecs')
    for ikr = 1:num_ir
        klist = vector_all[ikr]
        w = wk[ikr]
        m = umap[ikr]
        Nk = length(klist)
        @printf("ikr = %d, map = %d, w = %f, Nk = %d\n", ikr, m, w, Nk)
        for ik = 1:Nk
            kidx = klist[ik]
            kpt = b1*kidx[1]/meshk[1] + b2*kidx[2]/meshk[2] + b3*kidx[3]/meshk[3]
            @printf("map = %d ik = %4d: %8.5f %8.5f %8.5f  w = %8.5f Npw = %8d\n",
                    m, ik, kpt[1], kpt[2], kpt[3], w,
                    calc_Ngw(15.0, LatVecs, kpt))
        end
    end

end


function calc_Ngw( ecutwfc, LatVecs, kcart::Array{Float64,1} )

    ecutrho = 4.0*ecutwfc
    RecVecs = 2*pi*inv(LatVecs')

    LatVecsLen = zeros(3)
    LatVecsLen[1] = norm(LatVecs[:,1])
    LatVecsLen[2] = norm(LatVecs[:,2])
    LatVecsLen[3] = norm(LatVecs[:,3])

    Ns1 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[1]/pi ) + 1
    Ns2 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[2]/pi ) + 1
    Ns3 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[3]/pi ) + 1

    Ns1 = good_fft_order(Ns1)
    Ns2 = good_fft_order(Ns2)
    Ns3 = good_fft_order(Ns3)

    Ns = (Ns1,Ns2,Ns3)

    gvec = PWDFT.init_gvec( Ns, RecVecs, ecutrho )  # this is wasteful, need to be passed as input to the function
    Ng = gvec.Ng

    Gk2 = zeros(Float64,Ng)
    Gk = zeros(Float64,3)
    G = gvec.G
    
    for ig = 1:Ng
        Gk[:] = G[:,ig] .+ kcart[:]
        Gk2[ig] = Gk[1]^2 + Gk[2]^2 + Gk[3]^2
    end
    return length( findall( 0.5*Gk2 .<= ecutwfc ) )

end

function test_Al()
    # Atoms
    ALAT = 7.6525970200
    atoms = Atoms( xyz_string_frac=
        """
        1

        Al  0.0   0.0   0.0
        """, in_bohr=true,
        LatVecs=gen_lattice_fcc(ALAT) )

    kpoints_Ngw(atoms, [8,8,8])
end

test_Al()