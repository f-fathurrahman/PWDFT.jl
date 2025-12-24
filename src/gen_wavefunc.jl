# -----------------------------------------------------------------------------
# Generate "zeros" wavefunction
# -----------------------------------------------------------------------------

function zeros_BlochWavefunc( Ham::Hamiltonian )
    return zeros_BlochWavefunc( Ham.pw, Ham.electrons )
end

function zeros_BlochWavefunc( pw::PWGrid, electrons::Electrons )
    Nspin = electrons.Nspin_channel
    Nstates = electrons.Nstates
    Npol = 1
    if electrons.noncollinear
        Npol = 2
    end
    return zeros_BlochWavefunc( pw, Nstates, Nspin, Npol = Npol )
end

function zeros_BlochWavefunc( pw::PWGrid, Nstates::Int64, Nspin::Int64; Npol = 1)
    Nkpt = pw.gvecw.kpoints.Nkpt
    Ngw = pw.gvecw.Ngw
    Nkspin = Nspin*Nkpt
    psiks = BlochWavefunc(undef,Nkspin)
    for ispin in 1:Nspin, ik in 1:Nkpt
        ikspin = ik + (ispin-1)*Nkpt
        psiks[ikspin] = zeros(ComplexF64, Npol*Ngw[ik], Nstates)
    end
	return psiks
end



# -----------------------------------------------------------------------------
# Generate random wavefunction
# -----------------------------------------------------------------------------

function rand_BlochWavefunc( Ham::Hamiltonian )
    Npol = 1
    if Ham.electrons.noncollinear
        Npol = 2
    end
    # Overlap is not needed
    if !Ham.need_overlap
        return rand_BlochWavefunc(
            Ham.pw, Ham.electrons.Nstates,
            Ham.electrons.Nspin_channel,
            Npol = Npol
        )
    end

    # Overlap is needed
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Ngw = Ham.pw.gvecw.Ngw
    Nspin = Ham.electrons.Nspin_channel
    Nstates = Ham.electrons.Nstates
    Nkspin = Nspin*Nkpt
    psiks = BlochWavefunc(undef,Nkspin)
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        ikspin = ik + (ispin-1)*Nkpt
        psiks[ikspin] = rand_Wavefunc(Ngw[ik], Nstates, Npol = Npol)
        ortho_sqrt!(Ham, psiks[ikspin])
    end
    return psiks
end

function rand_BlochWavefunc( pw::PWGrid, electrons::Electrons )
    Nspin = electrons.Nspin_channel
    Nstates = electrons.Nstates
    Npol = 1
    if Ham.electron.noncollinear
        Npol = 2
    end
    return rand_BlochWavefunc( pw, Nstates, Nspin, Npol = Npol )
end

function rand_BlochWavefunc( pw::PWGrid, Nstates::Int64, Nspin::Int64; Npol = 1 )
    Nkpt = pw.gvecw.kpoints.Nkpt
    Ngw = pw.gvecw.Ngw
    Nkspin = Nspin*Nkpt
    psiks = BlochWavefunc(undef,Nkspin)
    for ispin in 1:Nspin, ik in 1:Nkpt
        ikspin = ik + (ispin-1)*Nkpt
        psiks[ikspin] = rand_Wavefunc(Ngw[ik], Nstates, Npol = Npol)
    end
    return psiks
end

function rand_Wavefunc( Nbasis, Nstates; Npol=1 )
    return ortho_sqrt( rand(ComplexF64, Npol*Nbasis, Nstates) )
end

#
# Adapted from pwscf
#
function rand_wfc( Ham::Hamiltonian )

    Nstates = Ham.electrons.Nstates
    idx_gw2g = Ham.pw.gvecw.idx_gw2g
    G = Ham.pw.gvec.G
    k = Ham.pw.gvecw.kpoints.k
    Nspin = Ham.electrons.Nspin_channel
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Ngw = Ham.pw.gvecw.Ngw
    Npol = 1
    if Ham.electrons.noncollinear
        Npol = 2
    end
    psiks = zeros_BlochWavefunc(Ham)
    for ispin in 1:Nspin, ik in 1:Nkpt
        ikspin = ik + (ispin - 1)*Nkpt
        # reshape return a reference not a copy
        psi = reshape(psiks[ikspin], Ngw[ik], Npol, Nstates)
        for ist in 1:Nstates, igk in 1:Ngw[ik]
            ig = idx_gw2g[ik][igk]
            rr = rand()
            arg = 2*pi*rand()
            num = rr*cos(arg) + im*rr*sin(arg)
            Gk1 = G[1,ig] + k[1,ik]
            Gk2 = G[2,ig] + k[2,ik]
            Gk3 = G[3,ig] + k[3,ik]
            Gw2 = Gk1^2 + Gk2^2 + Gk3^2
            denum = Gw2 + 1.0
            psi[igk,ipol,ist] = num/denum
        end
        Ham.ispin = ispin
        Ham.ik = ik
        ortho_sqrt!(Ham, psiks[ikspin])
        # psiks now should be orthonormal w.r.t to op_S, if applicable
        Hr = Hermitian(psiks[ikspin]' * op_H(Ham, psiks[ikspin]))
        evals, evecs = eigen(Hr)

        # Need to set this?
        Ham.electrons.ebands[:,ik] = evals

        psiks[ikspin] *= evecs # rotate
    end

    return psiks

end
 
#function atomic_wfc(Ham::Hamiltonian)
#    return Ham.
#end