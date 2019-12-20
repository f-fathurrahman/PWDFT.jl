function kernel_copy_to_fft_grid_gw2r!( ist::Int64, idx_gw2r_ik, psi, ctmp )
    igk = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x

    if igk <= length(idx_gw2r_ik)
        ip = idx_gw2r_ik[igk]
        ctmp[ip,ist] = psi[igk,ist]
    end
    
    return
end


function kernel_copy_to_fft_grid_gw2r_1state!( idx_gw2r_ik, psi, ctmp )
    igk = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x

    if igk <= length(idx_gw2r_ik)
        ip = idx_gw2r_ik[igk]
        ctmp[ip] = psi[igk]
    end
    
    return
end


function kernel_copy_from_fft_grid_gw2r!( ist::Int64, idx_gw2r_ik, ctmp, psi )
    igk = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x

    if igk <= length(idx_gw2r_ik)
        ip = idx_gw2r_ik[igk]
        psi[igk,ist] = ctmp[ip,ist]
    end
    
    return
end


function kernel_copy_from_fft_grid_gw2r_1state!( idx_gw2r_ik, ctmp, psi )
    igk = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x

    if igk <= length(idx_gw2r_ik)
        ip = idx_gw2r_ik[igk]
        psi[igk] = ctmp[ip]
    end
    
    return
end