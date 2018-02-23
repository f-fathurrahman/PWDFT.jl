function ortho_gram_schmidt( vin::Array{Complex128,2} )
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
