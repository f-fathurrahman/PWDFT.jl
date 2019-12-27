import PWDFT: calc_epsxc_PBE, calc_Vxc_PBE


# epsxc is always of type CuArray{Float64,1}
# Vxc is always of type CuArray{Float64,2}

# Rhoe can be spinpol or not
function calc_epsxc_PBE( xc_calc::XCCalculator, pw::CuPWGrid, Rhoe::CuArray{Float64,2} )
    Npoints = size(Rhoe,1)
    epsxc = CuArrays.zeros(Float64,Npoints)
    calc_epsxc_PBE!( xc_calc, pw, Rhoe, epsxc )
    return epsxc
end

function calc_Vxc_PBE( xc_calc::XCCalculator, pw::CuPWGrid, Rhoe::CuArray{Float64,2} )
    Npoints = size(Rhoe,1)
    Nspin = size(Rhoe,2)
    Vxc = CuArrays.zeros(Float64, Npoints, Nspin)
    if Nspin == 1
        calc_Vxc_PBE!( xc_calc, pw, Rhoe[:,1], Vxc )
        return Vxc        
    else
        calc_Vxc_PBE!( xc_calc, pw, Rhoe, Vxc )
        return Vxc
    end
end

#
# Inplace version
#


function kernel_calc_v_3_squared!( v, v2 )
    ip = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x
    
    Npoints = size(v2,1)  # Not of v !!!

    if ip <= Npoints
        v2[ip] = v[1,ip]*v[1,ip] + v[2,ip]*v[2,ip] + v[3,ip]*v[3,ip]
    end
    return
end


function kernel_epsxc_PBE!( Rhoe, gRhoe2, epsxc )

    ip = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x
    Npoints = size(Rhoe,1)
    
    if ip <= Npoints

        ρ = Rhoe[ip]
        ∇ρ2 = gRhoe2[ip]

        ss_x, _ = cu_XC_x_slater( ρ )
        gss_x, _, _ = cu_XC_x_pbe( ρ, ∇ρ2 )
        
        ss_c, _ = cu_XC_c_pw( ρ )
        gss_c, _, _ = cu_XC_c_pbe( ρ, ∇ρ2 )

        epsxc[ip] = ss_x + ss_c + ( gss_x + gss_c )/ρ

    end

    return
end



function calc_epsxc_PBE!( xc_calc::XCCalculator, pw::CuPWGrid, Rhoe::CuArray{Float64,1}, epsxc )

    Npoints = size(Rhoe,1)

    # calculate gRhoe2
    gRhoe = op_nabla( pw, Rhoe )
    gRhoe2 = CuArrays.zeros( Float64, Npoints )
    
    Nthreads = 256
    Nblocks = ceil( Int64, Npoints/Nthreads )
    @cuda threads=Nthreads blocks=Nblocks kernel_calc_v_3_squared!( gRhoe, gRhoe2 )

    @cuda threads=Nthreads blocks=Nblocks kernel_epsxc_PBE!( Rhoe, gRhoe2, epsxc )

    return
end



function kernel_epsxc_PBE_spin!( Rhoe, gRhoe2, gRhoe2_up, gRhoe2_dn, epsxc )

    ip = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x
    Npoints = size(Rhoe,1)
    
    if ip <= Npoints

        ρ_up = Rhoe[ip,1]
        ρ_dn = Rhoe[ip,2]
        ρ = ρ_up + ρ_dn
        ζ = (ρ_up - ρ_dn)/ρ

        ss_x = cu_XC_x_slater_spin_E( ρ, ζ )
        ss_c = cu_XC_c_pw_spin_E( ρ, ζ )

        gss_xup = cu_XC_x_pbe_E( 2*ρ_up, 4*gRhoe2_up[ip] )
        gss_xdn = cu_XC_x_pbe_E( 2*ρ_dn, 4*gRhoe2_dn[ip] )

        gss_x = 0.5*( gss_xup + gss_xdn )

        gss_c = cu_XC_c_pbe_spin_E( ρ, ζ, gRhoe2[ip] )

        epsxc[ip] = ( ss_x + ss_c ) + ( gss_x + gss_c )/ρ

    end

    return
end



function calc_epsxc_PBE!( xc_calc::XCCalculator, pw::CuPWGrid, Rhoe::CuArray{Float64,2}, epsxc )

    Npoints = size(Rhoe,1)
    Nspin = size(Rhoe,2)
    if Nspin == 1
        calc_epsxc_PBE!( xc_calc, pw, Rhoe[:,1], epsxc )
        return
    end

    Rhoe_total = Rhoe[:,1] + Rhoe[:,2]

    # calculate gRhoe2
    gRhoe_up = op_nabla( pw, Rhoe[:,1] )
    gRhoe_dn = op_nabla( pw, Rhoe[:,2] )
    gRhoe = op_nabla( pw, Rhoe_total )
    #
    gRhoe2_up = CuArrays.zeros( Float64, Npoints )
    gRhoe2_dn = CuArrays.zeros( Float64, Npoints )
    gRhoe2 = CuArrays.zeros( Float64, Npoints )

    Nthreads = 256
    Nblocks = ceil( Int64, Npoints/Nthreads )

    @cuda threads=Nthreads blocks=Nblocks kernel_calc_v_3_squared!( gRhoe, gRhoe2 )
    @cuda threads=Nthreads blocks=Nblocks kernel_calc_v_3_squared!( gRhoe_up, gRhoe2_up )
    @cuda threads=Nthreads blocks=Nblocks kernel_calc_v_3_squared!( gRhoe_dn, gRhoe2_dn )

    @cuda threads=Nthreads blocks=Nblocks kernel_epsxc_PBE_spin!( Rhoe, gRhoe2, gRhoe2_up, gRhoe2_dn, epsxc )

    return

end



function kernel_Vxc_PBE!( Rhoe, gRhoe, gRhoe2, h, Vxc )

    ip = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x
    Npoints = size(Rhoe,1)
    
    if ip <= Npoints
        
        ρ = Rhoe[ip]

        _, vx = cu_XC_x_slater( ρ )
        _, v1x, v2x = cu_XC_x_pbe( ρ, gRhoe2[ip] )
        
        _, vc = cu_XC_c_pw( ρ )
        _, v1c, v2c = cu_XC_c_pbe( ρ, gRhoe2[ip] )

        Vxc[ip] = vx + vc + v1x + v1c
        
        v2xc = v2x + v2c
        for i in 1:3
            h[i,ip] = v2xc*gRhoe[i,ip]
        end

    end

    return
end


function calc_Vxc_PBE!( xc_calc::XCCalculator, pw::CuPWGrid, Rhoe::CuArray{Float64,1}, Vxc )
    
    Npoints = size(Rhoe,1)

    # calculate gRhoe2
    gRhoe = op_nabla( pw, Rhoe )
    gRhoe2 = CuArrays.zeros( Float64, Npoints )
    
    Nthreads = 256
    Nblocks = ceil( Int64, Npoints/Nthreads )
    
    @cuda threads=Nthreads blocks=Nblocks kernel_calc_v_3_squared!( gRhoe, gRhoe2 )

    h = CuArrays.zeros(Float64,3,Npoints)

    @cuda threads=Nthreads blocks=Nblocks kernel_Vxc_PBE!( Rhoe, gRhoe, gRhoe2, h, Vxc )

    dh = op_nabla_dot(pw, h)

    Vxc[:] = Vxc[:] - dh[:]

    return
end


function kernel_Vxc_PBE_spin!( Rhoe, gRhoe_up, gRhoe_dn, gRhoe2, gRhoe2_up, gRhoe2_dn, h_up, h_dn, Vxc )
0
    ip = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x
    Npoints = size(Rhoe,1)

    if ip <= Npoints
        ρ_up = Rhoe[ip,1]
        ρ_dn = Rhoe[ip,2]
        ρ = ρ_up + ρ_dn
        ζ = (ρ_up - ρ_dn)/ρ

        _, vxup, vxdn = cu_XC_x_slater_spin( ρ, ζ )
        _, vcup, vcdn = cu_XC_c_pw_spin( ρ, ζ )

        _, v1xup, v2xup = cu_XC_x_pbe( 2*ρ_up, 4*gRhoe2_up[ip] )
        _, v1xdn, v2xdn = cu_XC_x_pbe( 2*ρ_dn, 4*gRhoe2_dn[ip] )

        v2xup = 2.0 * v2xup
        v2xdn = 2.0 * v2xdn

        _, v1cup, v1cdn, v2c = cu_XC_c_pbe_spin( ρ, ζ, gRhoe2[ip] )
        v2cup = v2c
        v2cdn = v2c
        v2cud = v2c

        Vxc[ip,1] = vxup + vcup + v1xup + v1cup
        Vxc[ip,2] = vxdn + vcdn + v1xdn + v1cdn

        for i in 1:3
           grup = gRhoe_up[i,ip]
           grdn = gRhoe_dn[i,ip]
           h_up[i,ip] = ( v2xup + v2cup ) * grup + v2cud * grdn
           h_dn[i,ip] = ( v2xdn + v2cdn ) * grdn + v2cud * grup
        end
    end

    return

end



function calc_Vxc_PBE!( xc_calc::XCCalculator, pw::CuPWGrid, Rhoe::CuArray{Float64,2}, Vxc::CuArray{Float64,2} )
    
    Npoints = size(Rhoe,1)
    Nspin = size(Rhoe,2)
    if Nspin == 1
        calc_Vxc_PBE!( xc_calc, pw, Rhoe[:,1], Vxc )
        return
    end

    Rhoe_total = Rhoe[:,1] + Rhoe[:,2]
    # calculate gRhoe2
    gRhoe_up = op_nabla( pw, Rhoe[:,1] )
    gRhoe_dn = op_nabla( pw, Rhoe[:,2] )
    gRhoe = op_nabla( pw, Rhoe_total )

    gRhoe2_up = CuArrays.zeros( Float64, Npoints )
    gRhoe2_dn = CuArrays.zeros( Float64, Npoints )
    gRhoe2 = CuArrays.zeros( Float64, Npoints )

    Nthreads = 256
    Nblocks = ceil( Int64, Npoints/Nthreads )

    @cuda threads=Nthreads blocks=Nblocks kernel_calc_v_3_squared!( gRhoe, gRhoe2 )
    @cuda threads=Nthreads blocks=Nblocks kernel_calc_v_3_squared!( gRhoe_up, gRhoe2_up )
    @cuda threads=Nthreads blocks=Nblocks kernel_calc_v_3_squared!( gRhoe_dn, gRhoe2_dn )

    h_up = CuArrays.zeros(Float64,3,Npoints)
    h_dn = CuArrays.zeros(Float64,3,Npoints)

    dh_up = CuArrays.zeros(Float64,Npoints)
    dh_dn = CuArrays.zeros(Float64,Npoints)

    @cuda threads=Nthreads blocks=Nblocks kernel_Vxc_PBE_spin!( Rhoe, gRhoe_up, gRhoe_dn, gRhoe2, gRhoe2_up, gRhoe2_dn, h_up, h_dn, Vxc )

    dh_up[:] = op_nabla_dot(pw, h_up)
    dh_dn[:] = op_nabla_dot(pw, h_dn)

    Vxc[:,1] = Vxc[:,1] - dh_up[:]
    Vxc[:,2] = Vxc[:,2] - dh_dn[:]        

    return
end
