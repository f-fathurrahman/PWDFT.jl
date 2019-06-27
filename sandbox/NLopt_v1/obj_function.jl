function obj_function!(
    Ham::Hamiltonian,
    psiks::BlochWavefunc;
    skip_ortho=false #,
#    evals::Array{Float64,2},
#    kT::Float64
)

    #Nspin = Ham.electrons.Nspin
    #Nelectrons = Ham.electrons.Nelectrons
    #wk = Ham.pw.gvecw.kpoints.wk
    #Ham.electrons.Focc, E_fermi = calc_Focc( Nelectrons, wk, kT, evals, Nspin )
    #Entropy = calc_entropy( wk, kT, evals, E_fermi, Nspin )

    if !skip_ortho
        for i = 1:length(psiks)
            ortho_sqrt!(psiks[i])
        end
    end

    Rhoe = calc_rhoe( Ham, psiks )
    update!( Ham, Rhoe )
    
    Ham.energies = calc_energies( Ham, psiks )
    #Ham.energies.mTS = Entropy

    return sum( Ham.energies )

end