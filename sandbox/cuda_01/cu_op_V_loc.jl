function kernel_copy_to_fft_grid_gw2r( idx_gw2r_ik, psi, ctmp )
    igk = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x

    if igk <= length(idx_gw2r_ik)
        ctmp[idx_gw2r_ik[igk]] = psi[igk]
    end
    
    return
end
