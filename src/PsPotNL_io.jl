import Base: println
function println( pspotNL::PsPotNL )
    @printf("\n")
    @printf("NbetaNL = %d", pspotNL.NbetaNL)
end