using CUDAnative
using CuArrays
using PWDFT

include("CuPWGrid.jl")

function kernel_op_K( G, idx_gw2g_k, k1, k2, k3, psi, Kpsi )

    igk = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x
    #strd = blockDim().x * gridDim().x

    if igk <= length(idx_gw2g_k)
        @inbounds ig = idx_gw2g_k[igk]
        @inbounds Gw2 = (G[1,ig] + k1)^2 + (G[2,ig] + k2)^2 + (G[3,ig] + k3)^2
        @inbounds Kpsi[igk] = psi[igk]*Gw2
    end

    return
end


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