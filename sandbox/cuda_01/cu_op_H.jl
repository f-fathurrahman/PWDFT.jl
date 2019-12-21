import PWDFT: op_H

function op_H( Ham::CuHamiltonian, psi )
    if Ham.pspotNL.NbetaNL > 0
        return op_K( Ham, psi ) + op_V_loc( Ham, psi ) + op_V_Ps_nloc( Ham, psi )
    else
        return op_K( Ham, psi ) + op_V_loc( Ham, psi )
    end
end

# These methods are not widely tested
import Base: *
*( Ham::CuHamiltonian, psi ) = op_H( Ham, psi )
*( psi, Ham::CuHamiltonian ) = adjoint( op_H( Ham, adjoint(psi) ) ) # two-times adjoint?