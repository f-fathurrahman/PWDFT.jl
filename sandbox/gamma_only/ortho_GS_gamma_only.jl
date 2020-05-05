function ortho_GS_gamma_only( vin::Array{ComplexF64,2} )
    Ncol = size(vin)[2]
    v = copy(vin)
    for ii = 1:Ncol
        zz = 2*dot( v[:,ii], v[:,ii] ) - v[1,ii]*v[1,ii]
        v[:,ii] = v[:,ii]/sqrt(zz)
        for jj = ii+1 : Ncol
            puv = prj_gamma_only( v[:,ii], v[:,jj] )
            v[:,jj] = v[:,jj] - puv*v[:,ii]
        end
    end
    return v
end

function ortho_GS_gamma_only!( v::Array{ComplexF64,2} )
    Ncol = size(v)[2]
    for ii = 1:Ncol
        zz = 2*dot( v[:,ii], v[:,ii] ) - conj(v[1,ii])*v[1,ii]
        v[:,ii] = v[:,ii]/sqrt(zz)
        for jj = ii+1 : Ncol
            puv = prj_gamma_only( v[:,ii], v[:,jj] )
            #println("puv = ", puv)
            v[:,jj] = v[:,jj] - puv*v[:,ii]
        end
        zz = 2*dot( v[:,ii], v[:,ii] ) - conj(v[1,ii])*v[1,ii]
        println("zz = ", zz)
    end
    return
end

function prj_gamma_only( u, v )
    num = 2*dot( u, v ) - conj(u[1])*v[1]
    denum = 2*dot( u, u ) - conj(u[1])*u[1]
    return  num/denum
end
