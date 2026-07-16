function exx_set_symm( sym_info, Ns )

    Nsyms = sym_info.Nsyms

    nr1, nr2, nr3 = Ns
    Npoints = nr1*nr2*nr3
    rir = zeros(Int64, Npoints, Nsyms)
    
    ftau = zeros(Int64, 3, Nsyms)
    s_scaled = zeros(Int64, 3, 3, Nsyms)

    scale_sym_ops!(sym_info, nr1, nr2, nr3, s_scaled, ftau)
    for isym in 1:Nsyms
        for k in 1:nr3, j in 1:nr2, i in 1:nr1
            ri, rj, rk = rotate_grid_point(
                s_scaled[:,:,isym], ftau[:,isym],
                i, j, k, nr1, nr2, nr3 
            )
            ir = i + (j-1)*nr1 + (k-1)*nr1*nr2
            rir[ir,isym] = ri + (rj-1)*nr1 + (rk-1)*nr1*nr2
        end
    end
    return rir
end