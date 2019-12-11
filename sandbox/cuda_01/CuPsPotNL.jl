struct CuPsPotNL
    NbetaNL::Int64
    prj2beta::Array{Int64,4}
    betaNL::CuArray{ComplexF64,3}
end

function CuPsPotNL( pspNL::PsPotNL )
    return CuPsPotNL( pspNL.NbetaNL, pspNL.prj2beta, CuArray(pspNL.betaNL) )
end
