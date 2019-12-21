struct CuPsPotNL
    NbetaNL::Int64
    prj2beta::Array{Int64,4}
    betaNL::Array{CuArray{ComplexF64,2},1}
end

function CuPsPotNL( pspNL::PsPotNL )
    Nkpt = length( pspNL.betaNL )
    betaNL = Array{CuArray{ComplexF64,2},1}(undef, Nkpt)
    for ik = 1:Nkpt
        betaNL[ik] = CuArray( pspNL.betaNL[ik] )
    end
    return CuPsPotNL( pspNL.NbetaNL, pspNL.prj2beta, betaNL )
end


import PWDFT: calc_betaNL_psi

function calc_betaNL_psi(
    ik::Int64,
    betaNL::Array{CuArray{ComplexF64,2},1},
    psi::CuArray{ComplexF64,2}
)

    Nstates = size(psi)[2]
    NbetaNL = size(betaNL[1],2)

    betaNL_psi = CuArrays.zeros( ComplexF64, Nstates, NbetaNL )
    betaNL_psi[:,:] = conj( psi' * betaNL[ik] )
    return betaNL_psi
end

function calc_betaNL_psi(
    ik::Int64,
    betaNL::Array{CuArray{ComplexF64,2},1},
    psi::Array{ComplexF64,1}
)
    NbetaNL = size(betaNL[1],2)

    betaNL_psi = CuArrays.zeros( ComplexF64, NbetaNL )
    betaNL_psi[:] = conj( psi' * betaNL[ik] )
    return betaNL_psi
end

