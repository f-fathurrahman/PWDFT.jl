"""
Calculate XC potential using PBE functional.
This function need `pw` as the first argument because it
is needed to calculate the gradient and various quantities.
This is fallback for spin-unpolarized system.
"""
function calc_Vxc_PBE( pw::PWGrid, Rhoe::Array{Float64,1} )
    Npoints = size(Rhoe)[1]

    # calculate gRhoe2
    gRhoe = op_nabla( pw, Rhoe ) # gRhoe = ∇⋅Rhoe
    gRhoe2 = zeros( Float64, Npoints )
    for ip = 1:Npoints
        gRhoe2[ip] = dot( gRhoe[:,ip], gRhoe[:,ip] )
    end

    Vxc = zeros( Float64, Npoints )
    Vgxc = zeros( Float64, Npoints )
    #
    ccall( (:calc_Vxc_PBE, LIBXC_SO_PATH), Nothing,
           (Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
           Npoints, Rhoe, gRhoe2, Vxc, Vgxc )

    h = zeros(Float64,3,Npoints)
    for ip = 1:Npoints
        h[1,ip] = Vgxc[ip] * gRhoe[1,ip]
        h[2,ip] = Vgxc[ip] * gRhoe[2,ip]
        h[3,ip] = Vgxc[ip] * gRhoe[3,ip]
    end

    # div ( vgrho * gRhoe )
    divh = op_nabla_dot( pw, h )

    for ip = 1:Npoints 
        Vxc[ip] = Vxc[ip] - 2.0*divh[ip]
        #Vxc[ip] = Vxc[ip] - divh[ip]
    end

    return Vxc
end

"""
Calculate XC energy per particle using PBE functional.
This function need `pw` as the first argument because it
is needed to calculate the gradient and various quantities.
This is fallback for spin-unpolarized system.
"""
function calc_epsxc_PBE( pw::PWGrid, Rhoe::Array{Float64,1} )
    Npoints = size(Rhoe)[1]
    epsxc = zeros( Float64, Npoints )

    # calculate gRhoe2
    gRhoe = op_nabla( pw, Rhoe )
    gRhoe2 = zeros( Float64, Npoints )
    for ip = 1:Npoints
        gRhoe2[ip] = dot( gRhoe[:,ip], gRhoe[:,ip] )
    end
    #
    ccall( (:calc_epsxc_PBE, LIBXC_SO_PATH), Nothing,
           (Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
           Npoints, Rhoe, gRhoe2, epsxc )
    #
    return epsxc
end

#
# Spin polarized versions
#

function calc_epsxc_PBE( pw::PWGrid, Rhoe::Array{Float64,2} )

    Nspin = size(Rhoe)[2]
    if Nspin == 1
        return calc_epsxc_PBE( pw, Rhoe[:,1] )
    end

    Npoints = size(Rhoe)[1]
    epsxc = zeros( Float64, Npoints )

    # calculate gRhoe2
    gRhoe_up = op_nabla( pw, Rhoe[:,1] )
    gRhoe_dn = op_nabla( pw, Rhoe[:,2] )
    gRhoe2 = zeros( Float64, 3, Npoints )
    for ip = 1:Npoints
        gRhoe2[1,ip] = dot( gRhoe_up[:,ip], gRhoe_up[:,ip] )
        gRhoe2[2,ip] = dot( gRhoe_up[:,ip], gRhoe_dn[:,ip] )
        gRhoe2[3,ip] = dot( gRhoe_dn[:,ip], gRhoe_dn[:,ip] )
    end

    Rhoe_tmp = zeros(2,Npoints)
    Rhoe_tmp[1,:] = Rhoe[:,1]
    Rhoe_tmp[2,:] = Rhoe[:,2]

    ccall( (:calc_epsxc_PBE_spinpol, LIBXC_SO_PATH), Nothing,
           (Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
           Npoints, Rhoe_tmp, gRhoe2, epsxc )

    return epsxc
end

# spin-polarized version
function calc_Vxc_PBE( pw::PWGrid, Rhoe::Array{Float64,2} )

    Nspin = size(Rhoe)[2]
    if Nspin == 1
        return calc_Vxc_PBE( pw, Rhoe[:,1] )
    end

    Npoints = size(Rhoe)[1]

    # calculate gRhoe2
    gRhoe_up = op_nabla( pw, Rhoe[:,1] ) # gRhoe = ∇⋅Rhoe
    gRhoe_dn = op_nabla( pw, Rhoe[:,2] )

    gRhoe2 = zeros( Float64, 3, Npoints )
    for ip = 1:Npoints
        gRhoe2[1,ip] = dot( gRhoe_up[:,ip], gRhoe_up[:,ip] )
        gRhoe2[2,ip] = dot( gRhoe_up[:,ip], gRhoe_dn[:,ip] )
        gRhoe2[3,ip] = dot( gRhoe_dn[:,ip], gRhoe_dn[:,ip] )
    end

    Vxc_tmp = zeros( Float64, 2*Npoints )
    Vgxc = zeros( Float64, 3, Npoints )
    Vxc = zeros( Float64, Npoints, 2 )
    
    Rhoe_tmp = zeros(2,Npoints)
    Rhoe_tmp[1,:] = Rhoe[:,1]
    Rhoe_tmp[2,:] = Rhoe[:,2]

    #
    ccall( (:calc_Vxc_PBE_spinpol, LIBXC_SO_PATH), Nothing,
           (Int64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
           Npoints, Rhoe_tmp, gRhoe2, Vxc_tmp, Vgxc )

    ipp = 0
    for ip = 1:2:2*Npoints
        ipp = ipp + 1
        Vxc[ipp,1] = Vxc_tmp[ip]
        Vxc[ipp,2] = Vxc_tmp[ip+1]
    end

    h = zeros(Float64,3,Npoints)

    #
    # spin up
    #
    for ip = 1:Npoints
        h[1,ip] = Vgxc[ip] * gRhoe_up[1,ip]
        h[2,ip] = Vgxc[ip] * gRhoe_up[2,ip]
        h[3,ip] = Vgxc[ip] * gRhoe_up[3,ip]
    end
    #
    divh = op_nabla_dot( pw, h )
    # spin up
    for ip = 1:Npoints 
        Vxc[ip,1] = Vxc[ip,1] - 2.0*divh[ip]
    end

    #
    # Spin down
    #
    for ip = 1:Npoints
        h[1,ip] = Vgxc[ip] * gRhoe_dn[1,ip]
        h[2,ip] = Vgxc[ip] * gRhoe_dn[2,ip]
        h[3,ip] = Vgxc[ip] * gRhoe_dn[3,ip]
    end
    #
    divh = op_nabla_dot( pw, h )
    # spin down
    for ip = 1:Npoints 
        Vxc[ip,2] = Vxc[ip,2] - 2.0*divh[ip]
    end

    return Vxc

end

