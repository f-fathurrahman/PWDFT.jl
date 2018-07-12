# psi can be either single or multicolumn
# This is handled by individual op_* functions
function op_H( Ham::Hamiltonian, psi )
    if Ham.pspotNL.NbetaNL > 0
        return op_K( Ham, psi ) + op_V_loc( Ham, psi ) + op_V_Ps_nloc( Ham, psi )
    else
        return op_K( Ham, psi ) + op_V_loc( Ham, psi )
    end
end
