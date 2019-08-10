function grad_obj_function!(
    Ham::Hamiltonian,
    psiks::BlochWavefunc,
    g::BlochWavefunc;
    skip_ortho=false
)

    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt

    if !skip_ortho
        for i = 1:length(psiks)
            ortho_sqrt!(psiks[i])
        end
    end

    Rhoe = calc_rhoe( Ham, psiks )
    update!( Ham, Rhoe )

    for ispin = 1:Nspin, ik = 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        ikspin = ik + (ispin-1)*Nkpt
        g[ikspin] = calc_grad( Ham, psiks[ikspin] )
    end

    return

end

function precond_grad!( Ham, g, Kg )
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    for ispin = 1:Nspin, ik = 1:Nkpt
        ikspin = ik + (ispin-1)*Nkpt
        Kg[ikspin] = Kprec( ik, Ham.pw, g[ikspin] )
    end
    return
end
