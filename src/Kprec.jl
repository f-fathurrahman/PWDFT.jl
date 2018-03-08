function Kprec( pw::PWGrid, psi::Array{Complex128,2} )

    Ngwx  = size(psi)[1]
    Nstates = size(psi)[2]
    idx_gw2g = pw.gvecw.idx_gw2g
    Gw2 = pw.gvec.G2[idx_gw2g]
    Kpsi  = zeros( Complex128, Ngwx, Nstates )

    for ist = 1:Nstates
        for ig = 1:Ngwx
            Kpsi[ig,ist] = psi[ig,ist] / ( 1.0 + Gw2[ig] )
        end
    end
    return Kpsi
end

# For comparison with non-preconditioned
function Kprec(psi)
    return psi
end
