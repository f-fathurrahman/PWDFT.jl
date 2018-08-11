function rand_Wavefunc( Nbasis, Nstates )
    return ortho_sqrt( rand(ComplexF64,Nbasis,Nstates) )
end

function zeros_BlochWavefunc( pw::PWGrid, electrons::Electrons )
    Nkpt = pw.gvecw.kpoints.Nkpt
    Nspin = electrons.Nspin
    Nkspin = Nspin*Nkpt
    Ngw = pw.gvecw.Ngw
    Nstates = electrons.Nstates

    psiks = Array{Array{ComplexF64,2},1}(undef,Nkspin)

    for ispin = 1:Nspin
    for ik = 1:Nkpt
        ikspin = ik + (ispin-1)*Nkpt
        psiks[ikspin] = zeros(ComplexF64,Ngw[ik],Nstates)
    end
    end
	return psiks
end

function rand_BlochWavefunc( pw::PWGrid, electrons::Electrons )

    Nkpt = pw.gvecw.kpoints.Nkpt
    Nspin = electrons.Nspin
    Nkspin = Nspin*Nkpt
    Ngw = pw.gvecw.Ngw
    Nstates = electrons.Nstates

    psiks = Array{Array{ComplexF64,2},1}(undef,Nkspin)

    for ispin = 1:Nspin
    for ik = 1:Nkpt
        ikspin = ik + (ispin-1)*Nkpt
        psiks[ikspin] = rand_Wavefunc(Ngw[ik],Nstates)
    end
    end

	return psiks

end


