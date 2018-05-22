"""
Calculate XC potential using VWN functional.
This is fallback for spin-unpolarized system.
"""
function calc_Vxc_VWN( Rhoe::Array{Float64,1} )
    Npoints = size(Rhoe)[1]
    Vxc = zeros( Float64, Npoints )
    #
    ccall( (:calc_Vxc_VWN, LIBXC_SO_PATH), Void,
           (Int64, Ptr{Float64}, Ptr{Float64}),
           Npoints, Rhoe, Vxc )
    #
    return Vxc
end

"""
Calculate XC energy per particle using VWN functional.
This is fallback for spin-unpolarized system.
"""
function calc_epsxc_VWN( Rhoe::Array{Float64,1} )
    Npoints = size(Rhoe)[1]
    epsxc = zeros( Float64, Npoints )
    #
    ccall( (:calc_epsxc_VWN, LIBXC_SO_PATH), Void,
           (Int64, Ptr{Float64}, Ptr{Float64}),
           Npoints, Rhoe, epsxc )
    #
    return epsxc
end


"""
Calculate XC potential using VWN functional.
This function works for both spin-polarized and spin-unpolarized system.
"""
function calc_Vxc_VWN( Rhoe::Array{Float64,2} )

    Nspin = size(Rhoe)[2]
    if Nspin == 1
        return calc_Vxc_VWN( Rhoe[:,1] )
    end

    Npoints = size(Rhoe)[1]
    Vxc = zeros( Float64, Npoints, 2 )
    Vxc_tmp = zeros( Float64, 2*Npoints )
    
    # This choice is transposed version of Rhoe
    Rhoe_tmp = zeros(2,Npoints)
    Rhoe_tmp[1,:] = Rhoe[:,1]
    Rhoe_tmp[2,:] = Rhoe[:,2]

    ccall( (:calc_Vxc_VWN_spinpol, LIBXC_SO_PATH), Void,
           (Int64, Ptr{Float64}, Ptr{Float64}),
           Npoints, Rhoe_tmp, Vxc_tmp )

    ipp = 0
    for ip = 1:2:2*Npoints
        ipp = ipp + 1
        Vxc[ipp,1] = Vxc_tmp[ip]
        Vxc[ipp,2] = Vxc_tmp[ip+1]
    end
    return Vxc
end

"""
Calculate XC energy per particle using VWN functional.
This function works for both spin-polarized and spin-unpolarized system.
"""
function calc_epsxc_VWN( Rhoe::Array{Float64,2} )

    Nspin = size(Rhoe)[2]
    if Nspin == 1
        return calc_epsxc_VWN( Rhoe[:,1] )
    end

    Npoints = size(Rhoe)[1]
    epsxc = zeros( Float64, Npoints )

    Rhoe_tmp = zeros(2,Npoints)
    Rhoe_tmp[1,:] = Rhoe[:,1]
    Rhoe_tmp[2,:] = Rhoe[:,2]

    ccall( (:calc_epsxc_VWN_spinpol, LIBXC_SO_PATH), Void,
           (Int64, Ptr{Float64}, Ptr{Float64}),
           Npoints, Rhoe_tmp, epsxc )
    #
    return epsxc
end

#
# VWN parameterization of the exchange correlation energy
# Adaptep from Arias
#
function excVWN( n::Array{Float64} )
    Npoints = size(n)[1]
    # Constants
    X1 =  0.75*(3.0/(2.0*pi))^(2.0/3.0)
    A  =  0.0310907
    x0 = -0.10498
    b  =  3.72744
    c  =  12.9352
    Q = sqrt(4*c-b*b)
    X0 = x0*x0 + b*x0 + c
    #
    out = zeros(Npoints)
    for ip = 1:Npoints
        rs = (4.0*pi/3.0*n[ip])^(-1.0/3.0)
        #@printf("%8d %18.10f %18.10f\n", ip, n[ip], rs)
        x = sqrt(rs)
        X = x*x + b*x + c
        out[ip] = -X1/rs + A*( log(x*x/X) + 2*b/Q*atan(Q/(2*x+b))
                  -(b*x0)/X0*( log((x-x0)*(x-x0)/X) + 2*(2*x0+b)/Q*atan(Q/(2*x+b)) ))
    end
    return out
end

#
# d/dn deriv of VWN parameterization of the exchange correlation energy
#
function excpVWN( n::Array{Float64} )
    Npoints = size(n)[1]
    # Constants
    X1 =  0.75*(3.0/(2.0*pi))^(2.0/3.0)
    A  =  0.0310907
    x0 = -0.10498
    b  =  3.72744
    c  = 12.9352
    Q  = sqrt(4.0*c - b*b)
    X0 = x0*x0 + b*x0 + c
    out = zeros(Npoints)
    for ip = 1:Npoints
        rs = (4.0*pi/3.0*n[ip])^(-1.0/3.0)
        x = sqrt(rs)
        X = x*x + b*x + c
        dx = 0.5/x
        out[ip] = dx*( 2.0*X1 / (rs.*x)
                  + A*( 2.0/x - (2*x+b)/X - 4*b/(Q*Q + (2*x+b)*(2*x+b))
                  -(b*x0)/X0*( 2.0/(x-x0) - (2*x+b)/X - 4*(2*x0+b) /
                  (Q*Q + (2*x+b) * (2*x+b)) ) ))
        out[ip] = (-rs./(3*n[ip]))*out[ip]
    end
    return out
end
