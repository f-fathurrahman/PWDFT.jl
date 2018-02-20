# psi can be either single or multicolumn
function op_H( Ham::PWHamiltonian, psi )
    return op_K( Ham, psi ) + op_V_loc( Ham, psi )
end
