"""
Simple preconditioner based on kinetic energy.

- `ik`: index of kpoints

- `pw`: an instance of `PWGrid`

- `psi`: wave function
"""
function Kprec( ik::Int64, pw::PWGrid, psi::Array{ComplexF64,2} )

    Ngw_ik  = size(psi)[1]
    Nstates = size(psi)[2]
    idx_gw2g = pw.gvecw.idx_gw2g[ik]
    G = pw.gvec.G
    k = pw.gvecw.kpoints.k[:,ik]

    Kpsi  = zeros( ComplexF64, Ngw_ik, Nstates )
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

"""
No instance of `pwgrid` is passed, no preconditioning is used.
"""
function Kprec(psi)
    return psi
end
