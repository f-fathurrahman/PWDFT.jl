import Base: print
function print( io::IO, pspotNL::PsPotNL )
    @printf(io, "\n")
    @printf(io, "NbetaNL = %d", pspotNL.NbetaNL)
end
print( pspotNL::PsPotNL ) = print( stdout, pspotNL )