function grad_obj_function!(
    Ham::Hamiltonian,
    psiks_::BlochWavefunc,
    g::BlochWavefunc,
    Haux_::Array{Matrix{ComplexF64},1},
    g_Haux::Array{Matrix{ComplexF64},1};
    kT::Float64=1e-3,    
    skip_ortho=false
)

    psiks = copy(psiks_)
    Haux = copy(Haux_)

    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin
    Nkspin = Nkpt*Nspin
    Nstates = Ham.electrons.Nstates
    Nelectrons = Ham.electrons.Nelectrons    
    wk = Ham.pw.gvecw.kpoints.wk

    evals = zeros( Float64, Nstates, Nkspin )
    U_Haux = copy(Haux)
    for i in 1:Nkspin
        evals[:,i], U_Haux[i] = eigen( Haux[i] )
    end

    Focc_old = copy(Ham.electrons.Focc)

    Ham.electrons.Focc, E_fermi = calc_Focc( Nelectrons, wk, kT, evals, Nspin )
    Entropy = calc_entropy( wk, kT, evals, E_fermi, Nspin )

    if !skip_ortho
        for i = 1:length(psiks)
            ortho_sqrt!(psiks[i])
        end
    end

    # rotate psiks
    #for i = 1:Nkspin
    #    psiks[i] = psiks[i]*U_Haux[i]
    #end

    Rhoe_old = copy( Ham.rhoe )
    Rhoe = calc_rhoe( Ham, psiks )
    update!( Ham, Rhoe )

    for ispin = 1:Nspin, ik = 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        i = ik + (ispin-1)*Nkpt
        g[i], _, _, g_Haux[i] = calc_grad_Haux( Ham, psiks[i], evals[:,i], kT )
    end

    update!( Ham, Rhoe_old )
    Ham.electrons.Focc = Focc_old

    return

end