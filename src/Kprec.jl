"""
Simple preconditioner based on kinetic energy.

- `ik`: index of kpoints

- `pw`: an instance of `PWGrid`

- `psi`: wave function
"""
function Kprec( ik::Int64, pw::PWGrid, psi )

    Ngw_ik  = size(psi)[1]
    Nstates = size(psi)[2]
    idx_gw2g = pw.gvecw.idx_gw2g[ik]
    G = pw.gvec.G
    k = pw.gvecw.kpoints.k[:,ik]

    Kpsi  = zeros( ComplexF64, Ngw_ik, Nstates )

    for ist = 1:Nstates
        for igk = 1:Ngw_ik
            ig = idx_gw2g[igk]
            Gw2 = (G[1,ig] + k[1])^2 + (G[2,ig] + k[2])^2 + (G[3,ig] + k[3])^2
            Kpsi[igk,ist] = psi[igk,ist] / ( 1.0 + Gw2 )
        end
    end
    return Kpsi
end

function Kprec!( ik::Int64, pw::PWGrid, psi, Kpsi )
    Ngw_ik  = size(psi)[1]
    Nstates = size(psi)[2]
    idx_gw2g = pw.gvecw.idx_gw2g[ik]
    G = pw.gvec.G
    @views k = pw.gvecw.kpoints.k[:,ik]
    for ist = 1:Nstates
        for igk = 1:Ngw_ik
            ig = idx_gw2g[igk]
            Gw2 = (G[1,ig] + k[1])^2 + (G[2,ig] + k[2])^2 + (G[3,ig] + k[3])^2
            Kpsi[igk,ist] = psi[igk,ist] / ( 1.0 + Gw2 )
        end
    end
    return
end


function Kprec_inplace!( ik::Int64, pw::PWGrid, psi )
    Ngw_ik  = size(psi)[1]
    Nstates = size(psi)[2]
    idx_gw2g = pw.gvecw.idx_gw2g[ik]
    G = pw.gvec.G
    @views k = pw.gvecw.kpoints.k[:,ik]
    for ist = 1:Nstates
        for igk = 1:Ngw_ik
            ig = idx_gw2g[igk]
            Gw2 = (G[1,ig] + k[1])^2 + (G[2,ig] + k[2])^2 + (G[3,ig] + k[3])^2
            psi[igk,ist] = psi[igk,ist] / ( 1.0 + Gw2 )
        end
    end
    return
end


"""
No instance of `pwgrid` is passed, no preconditioning is used.
"""
function Kprec(psi)
    return psi
end
