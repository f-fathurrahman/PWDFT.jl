# psi can be either Array{ComplexF64,1} (single columns)
# or Array{ComplexF64,2} (multicolumn)
# This is handled by individual op_* functions
function op_H( Ham::Hamiltonian, psi )
    # XXX This can be simplified, just call op_H! ?
    if isnothing(Ham.exx)
        if Ham.pspotNL.NbetaNL > 0
            return op_K( Ham, psi ) + op_V_loc( Ham, psi ) + op_V_Ps_nloc( Ham, psi ) + op_Vtau( Ham, psi )
        else
            return op_K( Ham, psi ) + op_V_loc( Ham, psi ) + op_Vtau( Ham, psi )
        end
    else
        if Ham.exx.is_active
            if Ham.pspotNL.NbetaNL > 0
                return op_K( Ham, psi ) + op_V_loc( Ham, psi ) + op_V_Ps_nloc( Ham, psi ) + op_Vtau( Ham, psi ) + op_Vexx(Ham, psi)
            else
                return op_K( Ham, psi ) + op_V_loc( Ham, psi ) + op_Vtau( Ham, psi ) + op_Vexx(Ham, psi)
            end
        else
            if Ham.pspotNL.NbetaNL > 0
                return op_K( Ham, psi ) + op_V_loc( Ham, psi ) + op_V_Ps_nloc( Ham, psi ) + op_Vtau( Ham, psi )
            else
                return op_K( Ham, psi ) + op_V_loc( Ham, psi ) + op_Vtau( Ham, psi )
            end 
        end
    end
end

function zero_out!( psiks::BlochWavefunc )
    N = length(psiks)
    for i in 1:N
        zero_out!(psiks[i])
    end
    return
end

function zero_out!( psi::Array{ComplexF64,2} )
    N = length(psi)
    for i in 1:N
        psi[i] = 0.0 + im*0.0
    end
    return
end

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
    op_Vexx!(Ham, psi, Hpsi)
    return
end

import Base: *
*( Ham::Hamiltonian, psi ) = op_H( Ham, psi )
*( psi::LinearAlgebra.Adjoint{Complex{Float64},Array{Complex{Float64},2}},
   Ham::Hamiltonian ) = adjoint( op_H( Ham, adjoint(psi) ) )