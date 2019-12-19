# epsxc is always of type Array{Float64,1}
# Vxc is always of type Array{Float64,2}

# Rhoe can be spinpol or not
function calc_epsxc_PBE( xc_calc::XCCalculator, pw::PWGrid, Rhoe )
    Npoints = size(Rhoe,1)
    epsxc = zeros(Float64,Npoints)
    calc_epsxc_PBE!( xc_calc, pw, Rhoe, epsxc )
    return epsxc
end

function calc_Vxc_PBE( xc_calc::XCCalculator, pw::PWGrid, Rhoe )
    Npoints = size(Rhoe,1)
    Nspin = size(Rhoe,2)
    Vxc = zeros(Float64, Npoints, Nspin)
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

function calc_epsxc_PBE!( xc_calc::XCCalculator, pw::PWGrid, Rhoe::Array{Float64,1}, epsxc )

    Npoints = size(Rhoe,1)

    # calculate gRhoe2
    gRhoe = op_nabla( pw, Rhoe )
    gRhoe2 = zeros( Float64, Npoints )
    for ip = 1:Npoints
        gRhoe2[ip] = gRhoe[1,ip]*gRhoe[1,ip] + gRhoe[2,ip]*gRhoe[2,ip] + gRhoe[3,ip]*gRhoe[3,ip]
    end

    for ip in 1:Npoints
        ρ = Rhoe[ip]
        ∇ρ2 = gRhoe2[ip]

        ss_x, _ = XC_x_slater( ρ )
        gss_x, _, _ = XC_x_pbe( ρ, ∇ρ2 )
        
        ss_c, _ = XC_c_pw( ρ )
        gss_c, _, _ = XC_c_pbe( ρ, ∇ρ2 )

        epsxc[ip] = ss_x + ss_c + ( gss_x + gss_c )/ρ
    end

    return
end


function calc_epsxc_PBE!( xc_calc::XCCalculator, pw::PWGrid, Rhoe::Array{Float64,2}, epsxc )

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
    gRhoe2_up = zeros( Float64, Npoints )
    gRhoe2_dn = zeros( Float64, Npoints )
    gRhoe2 = zeros( Float64, Npoints )

    # This is ugly but can save memory allocation a little bit
    for ip = 1:Npoints
        gRhoe2_up[ip] = gRhoe_up[1,ip]*gRhoe_up[1,ip] + gRhoe_up[2,ip]*gRhoe_up[2,ip] + gRhoe_up[3,ip]*gRhoe_up[3,ip]
        gRhoe2_dn[ip] = gRhoe_dn[1,ip]*gRhoe_dn[1,ip] + gRhoe_dn[2,ip]*gRhoe_dn[2,ip] + gRhoe_dn[3,ip]*gRhoe_dn[3,ip]
        gRhoe2[ip] = gRhoe[1,ip]*gRhoe[1,ip] + gRhoe[2,ip]*gRhoe[2,ip] + gRhoe[3,ip]*gRhoe[3,ip]
    end

    for ip in 1:Npoints

        ρ_up = Rhoe[ip,1]
        ρ_dn = Rhoe[ip,2]
        ρ = ρ_up + ρ_dn
        ζ = (ρ_up - ρ_dn)/ρ

        ss_x = XC_x_slater_spin_E( ρ, ζ )
        ss_c = XC_c_pw_spin_E( ρ, ζ )

        gss_xup = XC_x_pbe_E( 2*ρ_up, 4*gRhoe2_up[ip] )
        gss_xdn = XC_x_pbe_E( 2*ρ_dn, 4*gRhoe2_dn[ip] )

        gss_x = 0.5 * (gss_xup + gss_xdn)

        gss_c = XC_c_pbe_spin_E( ρ, ζ, gRhoe2[ip] )

        epsxc[ip] = ( ss_x + ss_c ) + ( gss_x + gss_c )/ρ

    end

    return

end


function calc_Vxc_PBE!( xc_calc::XCCalculator, pw::PWGrid, Rhoe::Array{Float64,1}, Vxc )
    
    Npoints = size(Rhoe,1)

    # calculate gRhoe2
    gRhoe = op_nabla( pw, Rhoe )
    gRhoe2 = zeros( Float64, Npoints )
    for ip = 1:Npoints
        gRhoe2[ip] = gRhoe[1,ip]*gRhoe[1,ip] + gRhoe[2,ip]*gRhoe[2,ip] + gRhoe[3,ip]*gRhoe[3,ip]
    end

    # h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
    h = zeros(Float64,3,Npoints)
    dh = zeros(Float64,Npoints)

    for ip in 1:Npoints

        ρ = Rhoe[ip]

        _, vx = XC_x_slater( ρ )
        _, v1x, v2x = XC_x_pbe( ρ, gRhoe2[ip] )
        
        _, vc = XC_c_pw( ρ )
        _, v1c, v2c = XC_c_pbe( ρ, gRhoe2[ip] )

        Vxc[ip] = vx + vc + v1x + v1c
        
        v2xc = v2x + v2c
        for i in 1:3
            h[i,ip] = v2xc*gRhoe[i,ip]
        end

    end

    dh[:] = op_nabla_dot(pw, h)

    for ip in 1:Npoints
        Vxc[ip] = Vxc[ip] - dh[ip]
    end

    return
end


function calc_Vxc_PBE!( xc_calc::XCCalculator, pw::PWGrid, Rhoe::Array{Float64,2}, Vxc::Array{Float64,2} )
    
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

    gRhoe2_up = zeros( Float64, Npoints )
    gRhoe2_dn = zeros( Float64, Npoints )
    gRhoe2 = zeros( Float64, Npoints )

    # This is ugly but can save memory allocation a little bit
    for ip = 1:Npoints
        gRhoe2_up[ip] = gRhoe_up[1,ip]*gRhoe_up[1,ip] + gRhoe_up[2,ip]*gRhoe_up[2,ip] + gRhoe_up[3,ip]*gRhoe_up[3,ip]
        gRhoe2_dn[ip] = gRhoe_dn[1,ip]*gRhoe_dn[1,ip] + gRhoe_dn[2,ip]*gRhoe_dn[2,ip] + gRhoe_dn[3,ip]*gRhoe_dn[3,ip]
        gRhoe2[ip] = gRhoe[1,ip]*gRhoe[1,ip] + gRhoe[2,ip]*gRhoe[2,ip] + gRhoe[3,ip]*gRhoe[3,ip]
    end

    # h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
    h_up = zeros(Float64,3,Npoints)
    h_dn = zeros(Float64,3,Npoints)
    #
    dh_up = zeros(Float64,Npoints)
    dh_dn = zeros(Float64,Npoints)

    for ip in 1:Npoints

        ρ_up = Rhoe[ip,1]
        ρ_dn = Rhoe[ip,2]
        ρ = ρ_up + ρ_dn
        ζ = (ρ_up - ρ_dn)/ρ

        _, vxup, vxdn = XC_x_slater_spin( ρ, ζ )
        _, vcup, vcdn = XC_c_pw_spin( ρ, ζ )

        _, v1xup, v2xup = XC_x_pbe( 2*ρ_up, 4*gRhoe2_up[ip] )
        _, v1xdn, v2xdn = XC_x_pbe( 2*ρ_dn, 4*gRhoe2_dn[ip] )

        v2xup = 2.0 * v2xup
        v2xdn = 2.0 * v2xdn

        _, v1cup, v1cdn, v2c = XC_c_pbe_spin( ρ, ζ, gRhoe2[ip] )
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

    dh_up[:] = op_nabla_dot(pw, h_up)
    dh_dn[:] = op_nabla_dot(pw, h_dn)

    for ip in 1:Npoints
        Vxc[ip,1] = Vxc[ip,1] - dh_up[ip]
        Vxc[ip,2] = Vxc[ip,2] - dh_dn[ip]        
    end

    return
end
