function ortho_GS_gamma( vin::Array{ComplexF64,2} )
    Ncol = size(vin)[2]
    v = copy(vin)
    for ii = 1:Ncol
        zz = 2*dot( v[:,ii], v[:,ii] ) - conj(v[1,ii])*v[1,ii]
        v[:,ii] = v[:,ii]/sqrt(zz)
        for jj = ii+1 : Ncol
            puv = prj_gamma( v[:,ii], v[:,jj] )
            v[:,jj] = v[:,jj] - puv*v[:,ii]
        end
    end
    return v
end

function ortho_GS_gamma!( v::Array{ComplexF64,2} )
    #println()
    #println("Entering ortho_GS_gamma:")
    
    Ncol = size(v)[2]
    
    for ii = 1:Ncol
        
        #zz = 2*dot( v[:,ii], v[:,ii] ) - conj(v[1,ii])*v[1,ii]
        c = dot(v[:,ii], v[:,ii])
        zz = c + conj(c) - conj(v[1,ii])*v[1,ii]

        v[:,ii] = v[:,ii]/sqrt(zz)
        
        #println("ii = ", ii, " zz = ", zz)
        
        for jj = ii+1 : Ncol
            puv = prj_gamma( v[:,ii], v[:,jj] )
            #println("puv = ", puv)
            v[:,jj] = v[:,jj] - puv*v[:,ii]
        end
    end
    return
end

function prj_gamma( u, v )
    c = dot(u,v)
    num = c + conj(c) - conj(u[1])*v[1]
    
    c = dot(u,u)
    denum = c + conj(c) - conj(u[1])*u[1]
    
    return num/denum

    #num = 2*dot( u, v ) - conj(u[1])*v[1]
    #denum = 2*dot( u, u ) - conj(u[1])*u[1]
    #return real(num/denum)
end
