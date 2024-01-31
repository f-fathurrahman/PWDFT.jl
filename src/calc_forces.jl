function calc_forces(
    Ham::Hamiltonian,
    psiks::BlochWavefunc,
)

    F_NN = calc_forces_NN( Ham.pw, Ham.atoms )

    F_Ps_loc = calc_forces_Ps_loc( Ham )

    F_Ps_nloc = calc_forces_Ps_nloc( Ham, psiks )
    symmetrize_vector!( Ham.pw.LatVecs, Ham.sym_info, F_Ps_nloc )
    # early return in handled in symmetrize_vector!
    
    F_nlcc = calc_forces_nlcc(Ham)
    F_scf_corr = calc_forces_scf_corr(Ham)

    return F_NN + F_Ps_loc + F_Ps_nloc + F_nlcc + F_scf_corr
end