import PWDFT: op_K

# Factor of 0.5 is included
function kernel_op_K!( ist, G, idx_gw2g_k, k1, k2, k3, psi, Kpsi )

    igk = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x

    if igk <= length(idx_gw2g_k)
        ig = idx_gw2g_k[igk]
        Gw2 = (G[1,ig] + k1)^2 + (G[2,ig] + k2)^2 + (G[3,ig] + k3)^2
        Kpsi[igk,ist] = 0.5*psi[igk,ist]*Gw2
    end

    return
end

function op_K( Ham::CuHamiltonian, psi::CuArray{ComplexF64,2} )

    # Get global index of current k-point index
    ik = Ham.ik

    Nstates = size(psi)[2]

    Ngwx = Ham.pw.gvecw.Ngwx
    Ngw = Ham.pw.gvecw.Ngw
    idx_gw2g_k = Ham.pw.gvecw.idx_gw2g[ik]
    G = Ham.pw.gvec.G
    k = Ham.pw.gvecw.kpoints.k[:,ik]

    Kpsi = CuArrays.zeros( ComplexF64, size(psi) )

    Nthreads = 256
    Nblocks = ceil( Int64, Ngwx/Nthreads )

    for ist = 1:Nstates
        @cuda threads=Nthreads blocks=Nblocks kernel_op_K!( ist, G, idx_gw2g_k, k[1], k[2], k[3], psi, Kpsi )
    end

    return Kpsi # factor of 0.5 is already included
end
