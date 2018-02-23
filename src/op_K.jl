# Apply kinetic operator to wave function in reciprocal space

function op_K( Ham::PWHamiltonian, psi::Array{Complex128,2} )
    #

    out = Array{Complex128}(size(psi))
    Nstates = size(psi)[2]

    pw = Ham.pw
    Ngwx = pw.gvecw.Ngwx
    G2   = pw.gvec.G2[pw.gvecw.idx_gw2r]

    for is = 1:Nstates
        for ig = 1:Ngwx
            out[ig,is] = psi[ig,is]*G2[ig]
        end
    end

    return 0.5*out # two minus signs -> positive
end


function op_K( Ham::PWHamiltonian, psi::Array{Complex128,1} )
    #
    out = Array{Complex128}(size(psi))
    pw = Ham.pw
    Ngwx = pw.gvecw.Ngwx
    G2   = pw.gvec.G2[pw.gvecw.idx_gw2r]
    for ig = 1:Ngwx
        out[ig] = psi[ig]*G2[ig]
    end
    return 0.5*out # two minus signs
end
