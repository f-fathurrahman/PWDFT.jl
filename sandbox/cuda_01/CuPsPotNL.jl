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
