import PWDFT: ortho_gram_schmidt, ortho_gram_schmidt!

function ortho_gram_schmidt( vin::CuArray{ComplexF64,2} )
    Ncol = size(vin)[2]
    v = copy(vin)
    for ii = 1:Ncol
        zz = dot( v[:,ii], v[:,ii] )
        v[:,ii] = v[:,ii]/sqrt(zz)
        for jj = ii+1 : Ncol
            puv = prj_GS( v[:,ii], v[:,jj] )
            v[:,jj] = v[:,jj] - puv*v[:,ii]
        end
    end
    return v
end

function ortho_gram_schmidt!( v::CuArray{ComplexF64,2} )
    Ncol = size(v)[2]
    for ii = 1:Ncol
        zz = dot( v[:,ii], v[:,ii] )
        v[:,ii] = v[:,ii]/sqrt(zz)
        for jj = ii+1 : Ncol
            puv = prj_GS( v[:,ii], v[:,jj] )
            v[:,jj] = v[:,jj] - puv*v[:,ii]
        end
    end
    return v
end

function prj_GS( u, v )
    return dot( u, v )/dot( u, u )
end
