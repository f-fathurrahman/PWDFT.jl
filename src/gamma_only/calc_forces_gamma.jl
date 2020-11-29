# skip symmetrize_forces
function calc_forces(
    Ham::HamiltonianGamma,
    psis::BlochWavefuncGamma,
)
    F_NN = calc_forces_NN( Ham.pw, Ham.atoms )
    F_Ps_loc = calc_forces_Ps_loc( Ham )
    F_Ps_nloc = calc_forces_Ps_nloc( Ham, psis ) # no symmetrization
    return F_NN + F_Ps_loc + F_Ps_nloc
end