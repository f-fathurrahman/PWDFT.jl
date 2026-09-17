# psi can be either Array{ComplexF64,1} (single columns)
# or Array{ComplexF64,2} (multicolumn)
# This is handled by individual op_* functions
function op_H( Ham::Hamiltonian, psi )
    Hpsi = similar(psi)
    op_H!(Ham, psi, Hpsi)
    return Hpsi
end

#XXX Need this?
import Base: fill!
function fill!( psiks::BlochWavefunc, x )
    N = length(psiks)
    for i in 1:N
        fill!(psiks[i], x)
    end
end

function op_H!( Ham::Hamiltonian, psi, Hpsi )
    fill!(Hpsi, 0.0 + im*0.0)
    op_K!( Ham, psi, Hpsi )
    op_V_loc!( Ham, psi, Hpsi )
    if Ham.pspotNL.NbetaNL > 0
        op_V_Ps_nloc!( Ham, psi, Hpsi )
    end
    op_Vtau!(Ham, psi, Hpsi)
    if !isnothing(Ham.exx)
        op_Vexx!(Ham, psi, Hpsi)
    end
    return
end

import Base: *
*( Ham::Hamiltonian, psi ) = op_H( Ham, psi )
*( psi::LinearAlgebra.Adjoint{Complex{Float64},Array{Complex{Float64},2}},
   Ham::Hamiltonian ) = adjoint( op_H( Ham, adjoint(psi) ) )