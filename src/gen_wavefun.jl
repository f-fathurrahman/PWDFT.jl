
# XXX Need to be updated for spin-polarized case
function gen_rand_wavefun( pw::PWGrid, electrons::Electrons; seed=1234 )

    Nkpt = pw.gvecw.kpoints.Nkpt
    Ngw = pw.gvecw.Ngw
    Nstates = electrons.Nstates

    psik = Array{Array{ComplexF64,2},1}(Nkpt)

    srand(seed)

    for ik = 1:Nkpt
        psik[ik] = rand(Ngw[ik],Nstates) + im*rand(Ngw[ik],Nstates)
        psik[ik] = ortho_gram_schmidt( psik[ik] )
    end

	return psik

end


