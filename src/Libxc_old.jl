using Libxc

# From Libxc.jl v0.1.3
mutable struct XCFuncType
end

#const LIBXC5 = "/home/efefer/mysoftwares/libxc-5.0.0/lib/libxc.so.9"
#const LIBXC5 = "/home/efefer/mysoftwares/libxc-4.3.4/lib/libxc.so.5"
const LIBXC5 = Libxc.libxc

function Libxc_xc_func_alloc()
    ccall(
        (:xc_func_alloc, LIBXC5), Ptr{XCFuncType},
        ()
    )
end

function Libxc_xc_func_init( p::Ptr{XCFuncType}, functional::Int, nspin::Int )
    ccall(
        (:xc_func_init, LIBXC5), Cint,
        (Ptr{XCFuncType}, Cint, Cint),
        p, functional, nspin
    )
end

function Libxc_xc_func_end( p::Ptr{XCFuncType} )
    ccall(
        (:xc_func_end, LIBXC5), Cvoid,
        (Ptr{XCFuncType},), p
    )
end

function Libxc_xc_func_free( p::Ptr{XCFuncType} )
    ccall(
        (:xc_func_free, LIBXC5), Cvoid,
        (Ptr{XCFuncType},), p
    )
end

function Libxc_xc_lda_vxc!(
    p::Ptr{XCFuncType},
    np::Int,
    rho::Array{Float64, 1},
    vrho::Array{Float64, 1}
)
    ccall( (:xc_lda_vxc, LIBXC5), Cvoid,
           (Ptr{XCFuncType}, Cint, Ptr{Float64}, Ptr{Float64}),
           p, np, rho, vrho)

end

function Libxc_xc_lda_exc!(
    p::Ptr{XCFuncType},
    np::Int,
    rho::Array{Float64,1},
    zk::Array{Float64,1}
)
    ccall( (:xc_lda_exc, LIBXC5), Cvoid,
           (Ptr{XCFuncType}, Cint, Ptr{Float64}, Ptr{Float64}),
           p, np, rho, zk
    )
end

function Libxc_xc_lda_exc_vxc!(
    p::Ptr{XCFuncType},
    np::Int,
    rho::Array{Float64,1},
    zk::Array{Float64,1},
    vrho::Array{Float64,1}
)
    ccall(
        (:xc_lda_exc_vxc, LIBXC5), Cvoid,
        (Ptr{XCFuncType}, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
        p, np, rho, zk, vrho
    )
end


function Libxc_xc_gga_exc!(
    p::Ptr{XCFuncType},
    np::Int,
    rho::Array{Float64,1},
    sigma::Array{Float64,1},
    zk::Array{Float64,1}
)
    ccall(
        (:xc_gga_exc, LIBXC5), Cvoid,
        (Ptr{XCFuncType}, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
        p, np, rho, sigma, zk
    )
end


function Libxc_xc_gga_exc_vxc!(
    p::Ptr{XCFuncType},
    np::Int,
    rho::Array{Float64,1},
    sigma::Array{Float64,1},
    zk::Array{Float64,1},
    vrho::Array{Float64,1},
    vsigma::Array{Float64,1}
)
    ccall(
        (:xc_gga_exc_vxc, LIBXC5), Cvoid,
        (Ptr{XCFuncType}, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
        p, np, rho, sigma, zk, vrho, vsigma
    )
end


function Libxc_xc_gga_vxc!(
    p::Ptr{XCFuncType},
    np::Int,
    rho::Array{Float64,1},
    sigma::Array{Float64,1},
    vrho::Array{Float64,1},
    vsigma::Array{Float64,1}
)
    ccall(
        (:xc_gga_vxc, LIBXC5), Cvoid,
        (Ptr{XCFuncType}, Cint, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
        p, np, rho, sigma, vrho, vsigma
    )
end


# MetaGGA
function Libxc_xc_mgga_exc!(
    p::Ptr{XCFuncType},
    np::Int,
    rho::Array{Float64,1},
    sigma::Array{Float64,1},
    lapl::Array{Float64,1},
    tau::Array{Float64,1},
    exc::Array{Float64,1}
)
    ccall(
        (:xc_mgga_exc, LIBXC5), Cvoid,
        ( Ptr{XCFuncType}, Cint,
          Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
          Ptr{Float64}
        ),
        p, np,
        rho, sigma, lapl, tau,
        exc
    )
end

function Libxc_xc_mgga_vxc!(
    p::Ptr{XCFuncType},
    np::Int,
    rho::Array{Float64,1},
    sigma::Array{Float64,1},
    lapl::Array{Float64,1},
    tau::Array{Float64,1},
    vrho::Array{Float64,1},
    vsigma::Array{Float64,1},
    vlapl::Array{Float64,1},
    vtau::Array{Float64,1}
)
    ccall(
        (:xc_mgga_vxc, LIBXC5), Cvoid,
        ( Ptr{XCFuncType}, Cint,
          Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
          Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}
        ),
        p, np,
        rho, sigma, lapl, tau,
        vrho, vsigma, vlapl, vtau
    )
end

function Libxc_xc_mgga_exc_vxc!(
    p::Ptr{XCFuncType},
    np::Int,
    rho::Array{Float64,1},
    sigma::Array{Float64,1},
    lapl::Array{Float64,1},
    tau::Array{Float64,1},
    exc::Array{Float64,1},
    vrho::Array{Float64,1},
    vsigma::Array{Float64,1},
    vlapl::Array{Float64,1},
    vtau::Array{Float64,1}
)
    ccall(
        (:xc_mgga_exc_vxc, LIBXC5), Cvoid,
        ( Ptr{XCFuncType}, Cint,
          Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
          Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}
        ),
        p, np,
        rho, sigma, lapl, tau,
        exc, vrho, vsigma, vlapl, vtau
    )
end