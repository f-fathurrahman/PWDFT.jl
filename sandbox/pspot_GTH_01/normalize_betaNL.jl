function normalize_betaNL(
    pw::PWGrid,
    betaNL::Array{ComplexF64,3},
    kpoints::KPoints
)

    Npoints = prod(pw.Ns)
    NbetaNL = size(betaNL)[2]

    Nkpt = kpoints.Nkpt
    Ngw = pw.gvecw.Ngw

    #
    # Check normalization in real space and reciprocal space
    #
    ctmp = zeros( ComplexF64, Npoints )
    for ik = 1:Nkpt
        #
        idx_gw2r = pw.gvecw.idx_gw2r[ik]
        #
        for ibeta = 1:NbetaNL
            norm_G = dot( betaNL[1:Ngw[ik],ibeta,ik], betaNL[1:Ngw[ik],ibeta,ik] )
            betaNL[1:Ngw[ik],ibeta,ik] = 1/sqrt(norm_G)*betaNL[1:Ngw[ik],ibeta,ik]         
        end
    end

end