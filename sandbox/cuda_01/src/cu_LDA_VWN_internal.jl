import PWDFT: calc_epsxc_VWN, calc_Vxc_VWN

# epsxc is always of type CuArray{Float64,1}
# Vxc is always of type CuArray{Float64,2}

# Rhoe can be spinpol or not
function calc_epsxc_VWN( xc_calc::XCCalculator, Rhoe::CuArray{Float64,2} )
    
    Npoints = size(Rhoe,1)
    epsxc = CuArrays.zeros(Float64,Npoints)
    
    calc_epsxc_VWN!( xc_calc, Rhoe, epsxc )
    
    return epsxc

end

function calc_Vxc_VWN( xc_calc::XCCalculator, Rhoe::CuArray{Float64,2} )
    
    Npoints = size(Rhoe,1)
    Nspin = size(Rhoe,2)
    
    Vxc = CuArrays.zeros(Float64, Npoints, Nspin)

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


function kernel_epsxc_VWN!( Rhoe, epsxc )
    ip = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x
    Npoints = size(Rhoe,1)
    if ip <= Npoints
        ss_x, _ = cu_XC_x_slater( Rhoe[ip] )
        ss_c, _ = cu_XC_c_vwn( Rhoe[ip] )
        epsxc[ip] = ss_x + ss_c
    end
    return
end

# Spin unpolarized version
function calc_epsxc_VWN!( xc_calc::XCCalculator, Rhoe::CuArray{Float64,1}, epsxc::CuArray{Float64,1} )
    Npoints = size(Rhoe,1)
    
    Nthreads = 256
    Nblocks = ceil(Int64, Npoints/Nthreads)
    @cuda threads=Nthreads blocks=Nblocks kernel_epsxc_VWN!( Rhoe, epsxc )
    
    return
end


function kernel_Vxc_VWN!( Rhoe, Vxc )
    ip = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x
    Npoints = size(Rhoe,1)
    if ip <= Npoints
        _, vx = cu_XC_x_slater( Rhoe[ip] )
        _, vc = cu_XC_c_vwn( Rhoe[ip] )
        Vxc[ip] = vx + vc
    end
    return
end


# Spin unpolarized version
function calc_Vxc_VWN!( xc_calc::XCCalculator, Rhoe::CuArray{Float64,1}, Vxc )
    Npoints = size(Rhoe,1)

    Nthreads = 256
    Nblocks = ceil(Int64, Npoints/Nthreads)
    @cuda threads=Nthreads blocks=Nblocks kernel_Vxc_VWN!( Rhoe, Vxc )

    return
end



function kernel_epsxc_VWN_spin!( Rhoe, epsxc )
    ip = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x
    
    Npoints = size(Rhoe,1)
    if ip <= Npoints
        ρ = Rhoe[ip,1] + Rhoe[ip,2]
        ζ = (Rhoe[ip,1] - Rhoe[ip,2])/ρ
        
        ss_x, _, _ = cu_XC_x_slater_spin( ρ, ζ )
        ss_c, _, _ = cu_XC_c_vwn_spin( ρ, ζ )
        
        epsxc[ip] = ss_x + ss_c
    end
    
    return
end


# Spin-polarized version
function calc_epsxc_VWN!( xc_calc::XCCalculator, Rhoe::CuArray{Float64,2}, epsxc::CuArray{Float64,1} )    
    Npoints = size(Rhoe,1)
    Nspin = size(Rhoe,2)
    if Nspin == 1
        calc_epsxc_VWN!( xc_calc, Rhoe[:,1], epsxc )
        return
    end

    Nthreads = 256
    Nblocks = ceil(Int64, Npoints/Nthreads)
    @cuda threads=Nthreads blocks=Nblocks kernel_epsxc_VWN_spin!( Rhoe, epsxc )

    return
end


function kernel_Vxc_VWN_spin!( Rhoe, Vxc )
    
    ip = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x
    
    Npoints = size(Rhoe,1)
    if ip <= Npoints
        ρ = Rhoe[ip,1] + Rhoe[ip,2]
        ζ = (Rhoe[ip,1] - Rhoe[ip,2])/ρ

        _, vxup, vxdn = cu_XC_x_slater_spin( ρ, ζ )
        _, vcup, vcdn = cu_XC_c_vwn_spin( ρ, ζ )

        Vxc[ip,1] = vxup + vcup
        Vxc[ip,2] = vxdn + vcdn
    end
    
    return
end



function calc_Vxc_VWN!( xc_calc::XCCalculator, Rhoe::CuArray{Float64,2}, Vxc::CuArray{Float64,2} )
    Npoints = size(Rhoe,1)
    Nspin = size(Rhoe,2)
    if Nspin == 1
        calc_Vxc_VWN!( xc_calc, Rhoe[:,1], Vxc )
        return
    end

    Nthreads = 256
    Nblocks = ceil(Int64, Npoints/Nthreads)
    @cuda threads=Nthreads blocks=Nblocks kernel_Vxc_VWN_spin!( Rhoe, Vxc )

    return
end
