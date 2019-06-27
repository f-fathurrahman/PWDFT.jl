function grad_obj_function!(
    Ham::Hamiltonian,
    psiks::BlochWavefunc,
    g::BlochWavefunc;
    rhoe_symm::Union{Nothing,RhoeSymmetrizer}=nothing,    
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
    if rhoe_symm != nothing
        #if Ham.sym_info.Nsyms > 1
        symmetrize_rhoe!( Ham, rhoe_symm, Rhoe )
    end
    update!( Ham, Rhoe )

    for ispin = 1:Nspin
        for ik = 1:Nkpt
            Ham.ispin = ispin
            Ham.ik = ik
            ikspin = ik + (ispin-1)*Nkpt
            g[ikspin] = calc_grad( Ham, psiks[ikspin] )
        end
    end

    return

end