include("PWDFT_cuda.jl")

function main()

    pw = PWGrid(20.0, gen_lattice_fcc(10.0))
    
    cu_gvec = CuGVectors( pw.gvec )
    cu_gvecw = CuGVectorsW( pw.gvecw )

    ik = 1
    Ngwk = cu_gvecw.Ngw[ik]
    Nstates = 4

    psi = CuArrays.rand( ComplexF64, Ngwk, Nstates )
    Kpsi = CuArrays.zeros( ComplexF64, Ngwk, Nstates )

    Nthreads = 256 # Nthreads per block
    Nblocks = ceil(Int64, Ngwk/Nthreads)
    println("Nblocks = ", Nblocks)

    k1 = cu_gvecw.kpoints.k[1,ik]
    k2 = cu_gvecw.kpoints.k[2,ik]
    k3 = cu_gvecw.kpoints.k[3,ik]

    G = cu_gvec.G
    idx_gw2g_k = cu_gvecw.idx_gw2g[ik]

    for ist in 1:Nstates
        @cuda threads=Nthreads blocks=Nblocks kernel_op_K!( ist, G, idx_gw2g_k, k1, k2, k3, psi, Kpsi )
    end

    println("Pass here")

end

main()