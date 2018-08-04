
# XXX Need to be updated for spin-polarized case
function gen_rand_wavefun( pw::PWGrid, electrons::Electrons; seed=1234 )

    Nkpt = pw.gvecw.kpoints.Nkpt
    Nspin = electrons.Nspin
    Nkspin = Nspin*Nkpt
    Ngw = pw.gvecw.Ngw
    Nstates = electrons.Nstates

    psiks = Array{Array{ComplexF64,2},1}(undef,Nkspin)

    Random.seed!(seed)

    for ispin = 1:Nspin
    for ik = 1:Nkpt
        ikspin = ik + (ispin-1)*Nkpt
        psiks[ikspin] = ortho_sqrt(rand(ComplexF64,Ngw[ik],Nstates))
    end
    end

	return psiks

end


