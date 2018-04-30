# Apply kinetic operator to wave function in reciprocal space

function op_K( Ham::PWHamiltonian, psi::Array{Complex128,2} )
    #
    ik = Ham.ik

    Nstates = size(psi)[2]

    pw = Ham.pw
    Ngwx = pw.gvecw.Ngwx
    Ngw = pw.gvecw.Ngw
    idx_gw2g = pw.gvecw.idx_gw2g[ik]
    k = pw.gvecw.kpoints.k[:,ik]

    Gw = zeros(3)
    out = Array{Complex128}(size(psi))

    for ist = 1:Nstates
        for igk = 1:Ngw[ik]
            ig = idx_gw2g[igk]
            Gw[:] = pw.gvec.G[:,ig] + k[:]
            Gw2 = Gw[1]^2 + Gw[2]^2 + Gw[3]^2
            out[igk,ist] = psi[igk,ist]*Gw2
        end
    end

    return 0.5*out # two minus signs -> positive
end


function op_K( Ham::PWHamiltonian, psi::Array{Complex128,1} )
    #
    ik = Ham.ik

    pw = Ham.pw
    Ngwx = pw.gvecw.Ngwx
    Ngw = pw.gvecw.Ngw
    idx_gw2g = pw.gvecw.idx_gw2g[ik]
    k = pw.gvecw.kpoints.k[:,ik]

    Gw = zeros(3)
    out = Array{Complex128}(size(psi))

    for igk = 1:Ngw[ik]
        ig = idx_gw2g[igk]
        Gw[:] = pw.gvec.G[:,ig] + k[:]
        Gw2 = Gw[1]^2 + Gw[2]^2 + Gw2[3]^2
        out[igk] = psi[igk]*Gw2
    end

    return 0.5*out # two minus signs -> positive
end



"""
function op_K( Ham::PWHamiltonian, psi::Array{Complex128,1} )
    #
    out = Array{Complex128}(size(psi))
    pw = Ham.pw
    Ngwx = pw.gvecw.Ngwx
    idx_gw2g = pw.gvecw.idx_gw2g
    Gw2   = pw.gvec.G2[idx_gw2g]

    for ig = 1:Ngwx
        out[ig] = psi[ig]*Gw2[ig]
    end

    return 0.5*out # two minus signs
end
"""