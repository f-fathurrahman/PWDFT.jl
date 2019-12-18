import PWDFT: calc_epsxc_VWN, calc_Vxc_VWN

# epsxc is always of type Array{Float64,1}
# Vxc is always of type Array{Float64,2}

# Rhoe can be spinpol or not
function calc_epsxc_VWN( xc_calc::XCCalculator, Rhoe )
    Npoints = size(Rhoe,1)
    epsxc = zeros(Float64,Npoints)
    calc_epsxc_VWN!( xc_calc, Rhoe, epsxc )
    return epsxc
end

function calc_Vxc_VWN( xc_calc::XCCalculator, Rhoe )
    Npoints = size(Rhoe,1)
    Nspin = size(Rhoe,2)
    Vxc = zeros(Float64, Npoints, Nspin)
    if Nspin == 1
        calc_Vxc_VWN!( xc_calc, Rhoe[:,1], Vxc )
        return Vxc        
    else
        calc_Vxc_VWN!( xc_calc, Rhoe, Vxc )
        return Vxc
    end

end



#
# Inplace version
#

# Spin unpolarized version
function calc_epsxc_VWN!( xc_calc::XCCalculator, Rhoe::Array{Float64,1}, epsxc )
    Npoints = size(Rhoe,1)    
    for ip in 1:Npoints
        ss_x, _ = XC_x_slater( Rhoe[ip] )
        ss_c, _ = XC_c_vwn( Rhoe[ip] )
        epsxc[ip] = ss_x + ss_c
    end
    return
end

# Spin unpolarized version
function calc_Vxc_VWN!( xc_calc::XCCalculator, Rhoe::Array{Float64,1}, Vxc )
    Npoints = size(Rhoe,1)
    for ip in 1:Npoints
        _, vx = XC_x_slater( Rhoe[ip] )
        _, vc = XC_c_vwn( Rhoe[ip] )
        Vxc[ip] = vx + vc
    end
    return
end


# Spin-polarized version
function calc_epsxc_VWN!( xc_calc::XCCalculator, Rhoe::Array{Float64,2}, epsxc::Array{Float64,1} )    
    Npoints = size(Rhoe,1)
    Nspin = size(Rhoe,2)
    if Nspin == 1
        calc_epsxc_VWN!( xc_calc, Rhoe[:,1], epsxc )
        return
    end

    for ip in 1:Npoints
        ρ = Rhoe[ip,1] + Rhoe[ip,2]
        ζ = (Rhoe[ip,1] - Rhoe[ip,2])/ρ
        
        ss_x, _, _ = XC_x_slater_spin( ρ, ζ )
        ss_c, _, _ = XC_c_vwn_spin( ρ, ζ )
        
        epsxc[ip] = ss_x + ss_c
    end
    return
end

function calc_Vxc_VWN!( xc_calc::XCCalculator, Rhoe::Array{Float64,2}, Vxc::Array{Float64,2} )
    Npoints = size(Rhoe,1)
    Nspin = size(Rhoe,2)
    if Nspin == 1
        calc_Vxc_VWN!( xc_calc, Rhoe[:,1], Vxc )
        return
    end
    for ip in 1:Npoints
        ρ = Rhoe[ip,1] + Rhoe[ip,2]
        ζ = (Rhoe[ip,1] - Rhoe[ip,2])/ρ

        _, vxup, vxdn = XC_x_slater_spin( ρ, ζ )
        _, vcup, vcdn = XC_c_vwn_spin( ρ, ζ )

        Vxc[ip,1] = vxup + vcup
        Vxc[ip,2] = vxdn + vcdn
    end    
    return
end