function Kprec( ik::Int64, pw::PWGrid, psi::Array{Complex128,2} )

    Ngw_ik  = size(psi)[1]
    Nstates = size(psi)[2]
    idx_gw2g = pw.gvecw.idx_gw2g[ik]
    G = pw.gvec.G
    k = pw.gvecw.kpoints.k[:,ik]

    Kpsi  = zeros( Complex128, Ngw_ik, Nstates )
    Gw = zeros(Float64,3)

    for ist = 1:Nstates
        for igk = 1:Ngw_ik
            ig = idx_gw2g[igk]
            Gw = G[:,ig] + k[:]
            Gw2 = Gw[1]^2 + Gw[2]^2 + Gw[3]^2
            Kpsi[igk,ist] = psi[igk,ist] / ( 1.0 + Gw2 )
        end
    end
    return Kpsi
end

# For comparison with non-preconditioned
function Kprec(psi)
    return psi
end
