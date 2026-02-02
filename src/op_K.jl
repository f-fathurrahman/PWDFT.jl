function op_K( Ham::Hamiltonian, psiks::BlochWavefunc )
    #    
    Nspin = Ham.electrons.Nspin_wf
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    out = zeros_BlochWavefunc(Ham)
    #
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
    Nspin = Ham.electrons.Nspin_wf
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
function op_K(
    Ham::Hamiltonian,
    psi::AbstractArray{ComplexF64}
)
    out = zeros(ComplexF64,size(psi))
    op_K!(Ham, psi, out)
    return out
end

# with preallocated array
# NOTE: The result is ACCUMULATED in Hpsi
function op_K!(
    Ham::Hamiltonian,
    psi::AbstractArray{ComplexF64},
    Hpsi::AbstractArray{ComplexF64}
)
    #
    Npol = 1
    if Ham.electrons.noncollinear
        Npol = 2
    end
    ik = Ham.ik # current ik index
    Nstates = size(psi,2)
    pw = Ham.pw
    Ngw = pw.gvecw.Ngw
    idx_gw2g = pw.gvecw.idx_gw2g[ik]
    G = pw.gvec.G
    k1 = pw.gvecw.kpoints.k[1,ik]
    k2 = pw.gvecw.kpoints.k[2,ik]
    k3 = pw.gvecw.kpoints.k[3,ik]
    #
    psir = reshape(psi, Ngw[ik], Npol, Nstates)
    Hpsir = reshape(Hpsi, Ngw[ik], Npol, Nstates)
    #
    for ist in 1:Nstates
        for ipol in 1:Npol, igk in 1:Ngw[ik]
            ig = idx_gw2g[igk]
            Gw2 = (G[1,ig] + k1)^2 + (G[2,ig] + k2)^2 + (G[3,ig] + k3)^2
            Hpsir[igk,ipol,ist] += 0.5*psir[igk,ipol,ist]*Gw2
        end
    end
    return
end
