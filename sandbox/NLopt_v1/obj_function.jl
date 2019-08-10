function obj_function!(
    Ham::Hamiltonian,
    psiks::BlochWavefunc;
    skip_ortho=false
)
    if !skip_ortho
        for i = 1:length(psiks)
            ortho_sqrt!(psiks[i])
        end
    end

    Rhoe = calc_rhoe( Ham, psiks )
    update!( Ham, Rhoe )
    
    Ham.energies = calc_energies( Ham, psiks )

    return sum( Ham.energies )

end