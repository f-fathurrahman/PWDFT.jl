function precKerker( pw::PWGrid, Nelectrons::Float64, R::Array{Float64,1} )
    Rg = R_to_G(pw, R)
    Ng = pw.gvec.Ng
    G2 = pw.gvec.G2
    rs = (3.0*pw.CellVolume/(4*pi)/Nelectrons)^(1.0/3.0)
    Ag = (12.0/pi)^(2.0/3.0)/rs
    idx_g2r = pw.gvec.idx_g2r
    for ig = 1:Ng
        ip = idx_g2r[ig]
        Rg[ip] = G2[ig]/(G2[ig] + Ag)*Rg[ip]
    end
    return real(G_to_R(pw, Rg))
end