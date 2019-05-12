function init_Gshells( gvec::GVectors )

    Ng = gvec.Ng

    eps8 = 1e-8

    G2_sorted = copy(gvec.G2)
    ngl = 1
    for ig = 2:Ng
        if G2_sorted[ig] > (G2_sorted[ig-1] + eps8)
            ngl = ngl + 1
        end
    end

    G2_gshells = zeros(ngl)
    idx_gshells = zeros(Int64,Ng)
    
    G2_gshells[1] = G2_sorted[1]
    idx_gshells[1] = 1

    igl = 1
    for ig = 2:Ng
        if G2_sorted[ig] > (G2_sorted[ig-1] + eps8)
            igl = igl + 1
            G2_gshells[igl] = G2_sorted[ig]
        end
        idx_gshells[ig] = igl
    end

    return G2_gshells, idx_gshells

end
