function precKerker( pw::PWGrid, R::Array{Float64,1} )
    Rg = R_to_G(pw, R)
    Ng = pw.gvec.Ng
    G2 = pw.gvec.G2
    idx_g2r = pw.gvec.idx_g2r
    for ig = 1:Ng
        ip = idx_g2r[ig]
        Rg[ip] = G2[ig]/(G2[ig] + 1.0)*Rg[ip]
    end
    return real(G_to_R(pw, Rg))
end