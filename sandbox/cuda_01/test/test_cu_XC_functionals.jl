using CuArrays
using CUDAnative

include("cu_XC_functionals_internal.jl")

function kernel_LDA_nospin!( Rhoe )
    idx = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x
    N = length(Rhoe)
    if idx <= N
        ex, vx = cu_XC_x_slater( Rhoe[idx] )
        ec, vc = cu_XC_c_pw( Rhoe[idx] )
        ec, vc = cu_XC_c_vwn( Rhoe[idx] )
    end
    return
end

function kernel_LDA_spin!( Rhoe, zeta )
    idx = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x
    N = length(Rhoe)
    if idx <= N
        ex, vxup, vxdw = cu_XC_x_slater_spin( Rhoe[idx], zeta[idx] )
        ec, vcup, vcdw = cu_XC_c_pw_spin( Rhoe[idx], zeta[idx] )
        ec, vcup, vcdw = cu_XC_c_vwn_spin( Rhoe[idx], zeta[idx] )
    end
    return
end


function kernel_GGA_nospin!( Rhoe, gRhoe2 )
    idx = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x
    N = length(Rhoe)
    if idx <= N
        sx, v1x, v2x = cu_XC_x_pbe( Rhoe[idx], gRhoe2[idx] )
        sc, v1c, v2c = cu_XC_c_pbe( Rhoe[idx], gRhoe2[idx] )        
    end
    return
end


function kernel_GGA_spin!( Rhoe, zeta, gRhoe2 )
    idx = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x
    N = length(Rhoe)
    if idx <= N
        sc, v1cup, v1cdw, v2c = cu_XC_c_pbe_spin( Rhoe[idx], zeta[idx], gRhoe2[idx] )
    end
    return
end

function main()
    
    Npoints = 1000
    
    Rhoe = abs.( CuArrays.rand(Float64, Npoints) )
    zeta = CuArrays.rand(Float64, Npoints)

    gRhoe2 = abs.( CuArrays.rand(Float64, Npoints) )

    Nthreads = 256
    Nblocks = ceil(Int64, Npoints/Nthreads)

    #@cuda threads=Nthreads blocks=Nblocks kernel_LDA_nospin!( Rhoe )
    
    #@cuda threads=Nthreads blocks=Nblocks kernel_LDA_spin!( Rhoe, zeta )

    #@cuda threads=Nthreads blocks=Nblocks kernel_GGA_nospin!( Rhoe, gRhoe2 )

    @cuda threads=Nthreads blocks=Nblocks kernel_GGA_spin!( Rhoe, zeta, gRhoe2 )

    println("Pass here")
end

main()
