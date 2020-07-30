# Verbose version

function ortho_gram_schmidt( vin::Array{ComplexF64,2} )
    Ncol = size(vin)[2]
    v = copy(vin)
    for ii = 1:Ncol
        zz = dot( v[:,ii], v[:,ii] )
        v[:,ii] = v[:,ii]/sqrt(zz)
        for jj = ii+1 : Ncol
            puv = prj( v[:,ii], v[:,jj] )
            v[:,jj] = v[:,jj] - puv*v[:,ii]
        end
    end
    return v
end

function ortho_gram_schmidt!( v::Array{ComplexF64,2} )
    println()
    println("Entering ortho_gram_schmidt:")

    Ncol = size(v)[2]
    for ii = 1:Ncol
        
        zz = dot( v[:,ii], v[:,ii] )
        v[:,ii] = v[:,ii]/sqrt(zz)

        println("ii = ", ii, " zz = ", zz)

        for jj = ii+1 : Ncol

            puv = prj( v[:,ii], v[:,jj] )
            println("puv = ", puv)

            v[:,jj] = v[:,jj] - puv*v[:,ii]
        end
    end
    return v
end

function ortho_gram_schmidt( Ncol, vin )
    v = copy(vin)
    for ii = 1:Ncol
        zz = dot( v[:,ii], v[:,ii] )
        v[:,ii] = v[:,ii]/sqrt(zz)
        for jj = ii+1 : Ncol
            puv = prj( v[:,ii], v[:,jj] )
            v[:,jj] = v[:,jj] - puv*v[:,ii]
        end
    end
    return v
end

function ortho_gram_schmidt!( Ncol, v )
    for ii = 1:Ncol
        zz = dot( v[:,ii], v[:,ii] )
        v[:,ii] = v[:,ii]/sqrt(zz)
        for jj = ii+1 : Ncol
            puv = prj( v[:,ii], v[:,jj] )
            v[:,jj] = v[:,jj] - puv*v[:,ii]
        end
    end
end

function prj( u, v )
    return dot( u, v )/dot( u, u )
end
