#=

(6 October 2018)

The purpose of the codes in this file is to investigate a new 
way to handle XC correlation functional.

Saving pointers of XCFuncType does not imply good performance.
Only slight improvement over current implementation is observed.

Probably, we only need to save x_func_id and c_func_id in order
to open possibility to use of all LibXC functionals.
=#


using Printf
using BenchmarkTools

using PWDFT

const LIBXC = PWDFT.LIBXC

struct XCFunctional
    x_func_id::Int64
    c_func_id::Int64
end

function XCFunctional( ;
  x_func_id=0, c_func_id=0, xc_name="LDA_VWN", Nspin=1
)
    
    if (x_func_id == 0) && ( c_func_id == 0 )
        if xc_name == "LDA_VWN"
            x_func_id = 1
            c_func_id = 7
        else
            error( @sprintf("xc_name is not valid and func_id is not set") )
        end
    end

    #ccall( (:xc_family_from_id,LIBXC), Cint, (Cint,), x_func_id )

    return XCFunctional( x_func_id, c_func_id )
end


function calc_epsxc( xc_func::XCFunctional, Rhoe::Array{Float64,1} )

    Npoints = size(Rhoe)[1]
    Nspin = 1
    eps_x = zeros(Float64,Npoints)
    eps_c = zeros(Float64,Npoints)

    x_func_id = xc_func.x_func_id
    c_func_id = xc_func.c_func_id

    ptr = ccall( (:xc_func_alloc, LIBXC), Ptr{XCFuncType}, () )

    #
    # exchange part 
    #
    ccall( (:xc_func_init, LIBXC), Cvoid,
           (Ptr{XCFuncType}, Cint, Cint),
            ptr, x_func_id, Nspin)
    
    ccall( (:xc_lda_exc, LIBXC), Cvoid,
           (Ptr{XCFuncType}, Cint, Ptr{Float64}, Ptr{Float64}),
           ptr, Npoints, Rhoe, eps_x )
    
    ccall( (:xc_func_end, LIBXC), Cvoid, (Ptr{XCFuncType},), ptr )

    #
    # correlation part 
    #
    ccall( (:xc_func_init, LIBXC), Cvoid,
            ( Ptr{XCFuncType}, Cint, Cint),
            ptr, c_func_id, Nspin)
    
    ccall( (:xc_lda_exc, LIBXC), Cvoid,
           (Ptr{XCFuncType}, Cint, Ptr{Float64}, Ptr{Float64}),
           ptr, Npoints, Rhoe, eps_c )
    
    ccall( (:xc_func_end, LIBXC), Cvoid, (Ptr{XCFuncType},), ptr )

    #
    ccall( (:xc_func_free, LIBXC), Cvoid, (Ref{XCFuncType},), ptr )

    return eps_x + eps_c

end


function test_small()
    Npoints = 5
    Rhoe = Array{Float64}(undef,Npoints)
    Rhoe[:] = [0.1, 0.2, 0.3, 0.4, 0.5]

    epsxc = calc_epsxc_VWN( Rhoe )
    Vxc = calc_Vxc_VWN( Rhoe )

    xc_func = XCFunctional()
    epsxc_v2 = calc_epsxc( xc_func, Rhoe )

    for ip = 1:Npoints
        @printf("%3d %18.10f %18.10f %18.10f\n", ip, Rhoe[ip], epsxc[ip], epsxc_v2[ip])
    end
end


function bench_small()
    Npoints = 5
    Rhoe = Array{Float64}(undef,Npoints)
    Rhoe[:] = [0.1, 0.2, 0.3, 0.4, 0.5]

    @btime calc_epsxc_VWN( $Rhoe )

    xc_func = XCFunctional()
    
    @btime calc_epsxc( $xc_func, $Rhoe )
end


function bench_large()
    Npoints = 5000
    Rhoe = abs.( rand(Npoints) )

    @btime calc_epsxc_VWN( $Rhoe )

    xc_func = XCFunctional()
    
    @btime calc_epsxc( $xc_func, $Rhoe )
end


function test_constructor()
    xc_func = XCFunctional()
    println(xc_func)
end

test_constructor()

test_small()

bench_small()

bench_large()

