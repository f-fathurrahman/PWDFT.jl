function obj_function!(
    Ham::Hamiltonian,
    psiks_::BlochWavefunc,
    Haux_::Array{Matrix{ComplexF64},1};
    kT::Float64=1e-3,
    skip_ortho=false
)

    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin
    Nkspin = Nkpt*Nspin
    Nstates = Ham.electrons.Nstates
    Nelectrons = Ham.electrons.Nelectrons    
    wk = Ham.pw.gvecw.kpoints.wk

    evals = zeros( Float64, Nstates, Nkspin )
    
    psiks = copy(psiks_)
    Haux = copy(Haux_)
    U_Haux = copy(Haux_)

    for i in 1:Nkspin
        evals[:,i], U_Haux[i] = eigen( Haux[i] )
        Haux[i] = diagm( 0 => evals[:,i] )
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
    
    Ham.energies = calc_energies( Ham, psiks )
    Ham.energies.mTS = Entropy

    # set Rhoe back
    update!( Ham, Rhoe_old )
    Ham.electrons.Focc = copy(Focc_old)

    return sum( Ham.energies )

end