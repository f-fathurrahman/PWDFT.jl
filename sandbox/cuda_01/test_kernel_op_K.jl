using CUDAnative
using CuArrays
using PWDFT

include("CuPWGrid.jl")

function main()

    pw = PWGrid(20.0, gen_lattice_fcc(10.0))
    
    cu_gvec = CuGVectors( pw.gvec )
    cu_gvecw = CuGVectorsW( pw.gvecw )

    ik = 1
    Ngwk = cu_gvecw.Ngw[ik]
    psi = CuArray( rand(ComplexF64,Ngwk) )
    Kpsi = CuArrays.fill( 0.0 + im*0.0, Ngwk )

    Nthreads = 256 # Nthreads per block
    Nblocks = ceil(Int64, Ngwk/Nthreads)
    println("Nblocks = ", Nblocks)

    k1 = cu_gvecw.kpoints.k[1,ik]
    k2 = cu_gvecw.kpoints.k[2,ik]
    k3 = cu_gvecw.kpoints.k[3,ik]

    @cuda threads=Nthreads blocks=Nblocks kernel_op_K( cu_gvec.G, cu_gvecw.idx_gw2g[ik], k1, k2, k3, psi, Kpsi )

    println("Pass here")

end

main()