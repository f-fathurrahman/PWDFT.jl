# -----------------------------------------------------------------------------
# Generate "zeros" wavefunction
# -----------------------------------------------------------------------------

function zeros_BlochWavefunc( Ham::Hamiltonian )
    return zeros_BlochWavefunc( Ham.pw, Ham.electrons )
end

function zeros_BlochWavefunc( pw::PWGrid, electrons::Electrons )
    Nspin = electrons.Nspin
    Nstates = electrons.Nstates
    return zeros_BlochWavefunc( pw, Nstates, Nspin )
end

function zeros_BlochWavefunc( pw::PWGrid, Nstates::Int64, Nspin::Int64)
    Nkpt = pw.gvecw.kpoints.Nkpt
    Ngw = pw.gvecw.Ngw

    Nkspin = Nspin*Nkpt

    psiks = Array{Array{ComplexF64,2},1}(undef,Nkspin)

    for ispin = 1:Nspin
    for ik = 1:Nkpt
        ikspin = ik + (ispin-1)*Nkpt
        psiks[ikspin] = zeros(ComplexF64,Ngw[ik],Nstates)
    end
    end
	return psiks
end



# -----------------------------------------------------------------------------
# Generate random wavefunction
# -----------------------------------------------------------------------------

function rand_BlochWavefunc( Ham::Hamiltonian )
    return rand_BlochWavefunc( Ham.pw, Ham.electrons )
end

function rand_BlochWavefunc( pw::PWGrid, electrons::Electrons )
    Nspin = electrons.Nspin
    Nstates = electrons.Nstates
    return rand_BlochWavefunc( pw, Nstates, Nspin )
end

function rand_BlochWavefunc( pw::PWGrid, Nstates::Int64, Nspin::Int64 )
    Nkpt = pw.gvecw.kpoints.Nkpt
    Ngw = pw.gvecw.Ngw

    Nkspin = Nspin*Nkpt

    psiks = BlochWavefunc(undef,Nkspin)

    for ispin = 1:Nspin
    for ik = 1:Nkpt
        ikspin = ik + (ispin-1)*Nkpt
        psiks[ikspin] = rand_Wavefunc(Ngw[ik],Nstates)
    end
    end

	return psiks
end

function rand_Wavefunc( Nbasis, Nstates )
    return ortho_sqrt( rand(ComplexF64,Nbasis,Nstates) )
end

