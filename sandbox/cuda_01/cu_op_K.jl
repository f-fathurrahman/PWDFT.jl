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