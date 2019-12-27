# not fully implemented yet

function ortho_sqrt( X::CuArrays{ComplexF64,2} )
    # sqrtm operation
    F = eigen( psiks[i]' * psiks[i] )
    retmat = (F.vectors * Diagonal(F.values)) * F.vectors'
    #Udagger = inv( psiks[i]' * psiks[i] )

    # inv operation is not yet implemented
end