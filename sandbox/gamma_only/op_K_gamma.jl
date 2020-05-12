import PWDFT: op_K

function op_K( Ham::HamiltonianGamma, psis::BlochWavefuncGamma )

    Nstates = size(psis.data[1],2)
    Nspin = Ham.electrons.Nspin
    out = zeros_BlochWavefuncGamma(Ham)
    
    for ispin in 1:Nspin
        #Ham.ispin = ispin #XXX This is not needed, spin dependence is carried out in psis
        out.data[ispin] = op_K( Ham, psis.data[ispin] )
    end
    return out
end


# Apply kinetic operator to wave function in reciprocal space

function op_K( Ham::HamiltonianGamma, psi::Array{ComplexF64,2} )    
    Nstates = size(psi,2)
    Ngw = Ham.pw.gvecw.Ngw
    idx_gw2g = Ham.pw.gvecw.idx_gw2g
    G2 = Ham.pw.gvec.G2

    out = zeros(ComplexF64,size(psi))
    for ist in 1:Nstates, igw in 1:Ngw
        ig = idx_gw2g[igw]
        out[igw,ist] = psi[igw,ist]*G2[ig]
    end
    return 0.5*out # two minus signs -> positive
end


# This function is used by CheFSI
# FIXME: use multiple dispatch
function op_K( Ham::HamiltonianGamma, psi::Array{ComplexF64,1} )
    Nstates = size(psi,2)
    Ngw = Ham.pw.gvecw.Ngw
    idx_gw2g = Ham.pw.gvecw.idx_gw2g
    G2 = Ham.pw.gvec.G2
    out = zeros(ComplexF64,size(psi))
    for igw in 1:Ngw
        ig = idx_gw2g[igw]
        out[igw] = psi[igw]*G2[ig]
    end
    return 0.5*out # two minus signs -> positive
end
