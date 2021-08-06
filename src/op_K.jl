function op_K( Ham::Hamiltonian, psiks::BlochWavefunc )
    Nstates = size(psiks[1])[2] # Nstates should be similar for all Bloch states
    
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    out = zeros_BlochWavefunc(Ham)
    
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        out[ikspin] = op_K( Ham, psiks[ikspin] )
    end
    return out
end

# In-place, accumulated version
function op_K!( Ham::Hamiltonian, psiks::BlochWavefunc, Hpsiks::BlochWavefunc )
    #
    Nstates = size(psiks[1],2) # Nstates should be similar for all Bloch states
    #
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    #
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        i = ik + (ispin - 1)*Nkpt
        op_K!( Ham, psiks[i], Hpsiks[i] )
    end
    return
end


# Apply kinetic operator to wave function in reciprocal space

function op_K( Ham::Hamiltonian, psi::AbstractArray{ComplexF64,2} )
    #
    ik = Ham.ik

    Nstates = size(psi,2)

    pw = Ham.pw
    Ngwx = pw.gvecw.Ngwx
    Ngw = pw.gvecw.Ngw
    idx_gw2g = pw.gvecw.idx_gw2g[ik]
    G = pw.gvec.G
    k = pw.gvecw.kpoints.k[:,ik]

    out = zeros(ComplexF64,size(psi))

    for ist in 1:Nstates
        for igk in 1:Ngw[ik]
            ig = idx_gw2g[igk]
            Gw2 = (G[1,ig] + k[1])^2 + (G[2,ig] + k[2])^2 + (G[3,ig] + k[3])^2
            out[igk,ist] = psi[igk,ist]*Gw2
        end
    end

    return 0.5*out # two minus signs -> positive
end

# with preallocated array
# NOTE: The result is ACCUMULATED in Hpsi
function op_K!(
    Ham::Hamiltonian, psi::Array{ComplexF64,2},
    Hpsi::Array{ComplexF64,2}
)
    #
    ik = Ham.ik

    Nstates = size(psi)[2]

    pw = Ham.pw
    Ngwx = pw.gvecw.Ngwx
    Ngw = pw.gvecw.Ngw
    idx_gw2g = pw.gvecw.idx_gw2g[ik]
    G = pw.gvec.G
    k1 = pw.gvecw.kpoints.k[1,ik]
    k2 = pw.gvecw.kpoints.k[2,ik]
    k3 = pw.gvecw.kpoints.k[3,ik]
    for ist in 1:Nstates
        for igk in 1:Ngw[ik]
            ig = idx_gw2g[igk]
            Gw2 = (G[1,ig] + k1)^2 + (G[2,ig] + k2)^2 + (G[3,ig] + k3)^2
            Hpsi[igk,ist] = Hpsi[igk,ist] + 0.5*psi[igk,ist]*Gw2
        end
    end
    return
end



# This function is used by CheFSI
function op_K( Ham::Hamiltonian, psi::Array{ComplexF64,1} )
    #
    ik = Ham.ik

    pw = Ham.pw
    Ngwx = pw.gvecw.Ngwx
    Ngw = pw.gvecw.Ngw
    idx_gw2g = pw.gvecw.idx_gw2g[ik]
    G = pw.gvec.G
    k = pw.gvecw.kpoints.k[:,ik]

    out = zeros(ComplexF64,size(psi))

    for igk in 1:Ngw[ik]
        ig = idx_gw2g[igk]
        Gw2 = (G[1,ig] + k[1])^2 + (G[2,ig] + k[2])^2 + (G[3,ig] + k[3])^2
        out[igk] = psi[igk]*Gw2
    end

    return 0.5*out # two minus signs -> positive
end
