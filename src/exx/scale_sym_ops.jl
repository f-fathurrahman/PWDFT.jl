
function scale_sym_ops!(sym_info, nr1, nr2, nr3, s_scaled, ftau)
    # output: s_scaled and f_tau
    
    #=
    Generate rotation matrices and fractional translations in scaled
    crystallographic axis from symmetries, check their consistency.
    For use in real-space symmetrization routine "rotate_grid_point"
    =#

    Nsyms = sym_info.Nsyms
    s = sym_info.s
    ft = sym_info.ft

    SMALL = 1e-5
    ft_ = zeros(Float64, 3)
    n_bad = 0
    for isym in 1:Nsyms
        # check if rotation sends the FFT grid into itself
        if mod( s[2,1,isym]*nr1, nr2) != 0 ||
           mod( s[3,1,isym]*nr1, nr3) != 0 ||
           mod( s[1,2,isym]*nr2, nr1) != 0 ||
           mod( s[3,2,isym]*nr2, nr3) != 0 ||
           mod( s[1,3,isym]*nr3, nr1) != 0 ||
           mod( s[2,3,isym]*nr3, nr2) != 0
            warn("found rotation not compatible with FFT grid")
            bad = bad + 1
        end   
        s_scaled[1,1,isym] = s[1,1,isym]
        s_scaled[2,1,isym] = s[2,1,isym] * nr1 / nr2
        s_scaled[3,1,isym] = s[3,1,isym] * nr1 / nr3
        s_scaled[1,2,isym] = s[1,2,isym] * nr2 / nr1
        s_scaled[2,2,isym] = s[2,2,isym]
        s_scaled[3,2,isym] = s[3,2,isym] * nr2 / nr3
        s_scaled[1,3,isym] = s[1,3,isym] * nr3 / nr1
        s_scaled[2,3,isym] = s[2,3,isym] * nr3 / nr2
        s_scaled[3,3,isym] = s[3,3,isym]
        ft_[1] = ft[1,isym] * nr1
        ft_[2] = ft[2,isym] * nr2
        ft_[3] = ft[3,isym] * nr3
        #
        # check if the fractional translations are commensurate
        # with the FFT grid
        if any( abs.(ft_ - round.(Int64, ft_)) .> SMALL )
            warn("Found fractional translation not compatible with FFT grid")
            n_bad += 1
        end
        ftau[:,isym] = round.(Int64, ft_)
    end
    if n_bad > 0
        error("incompatible FFT grid")
    end
    return
end
