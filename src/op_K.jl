# Apply kinetic operator to wave function in reciprocal space

function op_K( Ham::PWHamiltonian, psi::Array{Complex128,2} )
    #

    out = Array{Complex128}(size(psi))
    Nstates = size(psi)[2]

    pw = Ham.pw
    Ngwx = pw.gvecw.Ngwx
    idx_gw2g = pw.gvecw.idx_gw2g
    Gw2   = pw.gvec.G2[idx_gw2g]

    for is = 1:Nstates
        for ig = 1:Ngwx
            out[ig,is] = psi[ig,is]*Gw2[ig]
        end
    end

    return 0.5*out # two minus signs -> positive
end


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
