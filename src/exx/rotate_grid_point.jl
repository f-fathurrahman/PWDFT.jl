function rotate_grid_point(s_scaled, ftau, i, j, k, nr1, nr2, nr3)
    #=
    This routine finds the symmetry-rotated point (ri, rj, rk)
    from a point (i, j, k) of the FFT real-space grid.
    The symmetry operation can include a fractional translation
    =#
    ri = s_scaled[1,1]*(i-1) + s_scaled[2,1]*(j-1) + s_scaled[3,1]*(k-1) - ftau[1]
    ri = mod(ri, nr1) + 1
    if ri < 1
        ri = ri + nr1
    end
    
    rj = s_scaled[1,2]*(i-1) + s_scaled[2,2]*(j-1) + s_scaled[3,2]*(k-1) - ftau[2]
    rj = mod(rj, nr2) + 1
    if rj < 1
        rj = rj + nr2
    end
    
    rk = s_scaled[1,3]*(i-1) + s_scaled[2,3]*(j-1) + s_scaled[3,3]*(k-1) - ftau[3]
    rk = mod(rk, nr3) + 1
    if rk < 1
        rk = rk + nr3
    end

    return ri, rj, rk
end
