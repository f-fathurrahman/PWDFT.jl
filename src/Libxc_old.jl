using Libxc

# From Libxc.jl v0.1.3
mutable struct XCFuncType
end

function Libxc_xc_func_alloc()
    ccall(
        (:xc_func_alloc, Libxc.libxc), Ptr{XCFuncType},
        ()
    )
end

function Libxc_xc_func_init( p::Ptr{XCFuncType}, functional::Int, nspin::Int )
    ccall(
        (:xc_func_init, Libxc.libxc), Cint,
        (Ptr{XCFuncType}, Cint, Cint),
        p, functional, nspin
    )
end

function Libxc_xc_func_end( p::Ptr{XCFuncType} )
    ccall(
        (:xc_func_end, Libxc.libxc), Cvoid,
        (Ptr{XCFuncType},), p
    )
end

function Libxc_xc_func_free( p::Ptr{XCFuncType} )
    ccall(
        (:xc_func_free, Libxc.libxc), Cvoid,
        (Ptr{XCFuncType},), p
    )
end

function Libxc_xc_lda_vxc!(
    p::Ptr{XCFuncType},
    np::Int,
    rho::Array{Float64, 1},
    vrho::Array{Float64, 1}
)
    ccall( (:xc_lda_vxc, Libxc.libxc), Cvoid,
           (Ptr{XCFuncType}, Cint, Ptr{Float64}, Ptr{Float64}),
           p, np, rho, vrho)

end

function Libxc_xc_lda_exc!(
    p::Ptr{XCFuncType},
    np::Int,
    rho::Array{Float64,1},
    zk::Array{Float64,1}
)
    ccall( (:xc_lda_exc, Libxc.libxc), Cvoid,
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
        (:xc_lda_exc_vxc, Libxc.libxc), Cvoid,
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
        (:xc_gga_exc, Libxc.libxc), Cvoid,
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
        (:xc_gga_exc_vxc, Libxc.libxc), Cvoid,
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
        (:xc_gga_vxc, Libxc.libxc), Cvoid,
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
        (:xc_mgga_exc, Libxc.libxc), Cvoid,
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
        (:xc_mgga_vxc, Libxc.libxc), Cvoid,
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
        (:xc_mgga_exc_vxc, Libxc.libxc), Cvoid,
        ( Ptr{XCFuncType}, Cint,
          Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
          Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}
        ),
        p, np,
        rho, sigma, lapl, tau,
        exc, vrho, vsigma, vlapl, vtau
    )
end