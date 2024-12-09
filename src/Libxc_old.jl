using Libxc_jll


# This is a hack for testing various difference between Libxc implementation
# I mainly used this for testing metaGGA functionals.
#
#const LIBXC5 = "/home/efefer/mysoftwares/libxc-5.0.0/lib/libxc.so.9"
#const LIBXC5 = "/home/efefer/mysoftwares/libxc-4.3.4dyn/lib/libxc.so.5"
#const LIBXC5 = "/usr/lib/x86_64-linux-gnu/libxc.so.5"
const LIBXC5 = Libxc_jll.libxc


# From Libxc.jl v0.1.3
mutable struct XCFuncType
end

mutable struct FuncReferType
    ref::Ptr{UInt8}
    doi::Ptr{UInt8}
    bibtex::Ptr{UInt8}
end

mutable struct FuncParamsType
    value::Float64
    description::Ptr{UInt8}
end

# FIXME: Need this?
# MAX_REFERENCES=5
mutable struct XCFuncInfoType
    number::Cint
    kind::Cint
    name::Ptr{UInt8}
    family::Cint
    refs::NTuple{5, Ref{FuncReferType}}
    #
    flags::Cint
    #
    dens_threshold::Float64
    #
    n_ext_params::Cint
    ext_params::Ref{FuncParamsType}
    set_ext_params::Ptr{Cvoid}
    #
    init::Ptr{Cvoid}
    end_::Ptr{Cvoid}
    lda::Ptr{Cvoid}
    gga::Ptr{Cvoid}
    mgga::Ptr{Cvoid}
end




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


function Libxc_xc_func_set_dens_threshold(p::Ptr{XCFuncType}, dens_threshold::Float64)
    ccall(
        (:xc_func_set_dens_threshold, LIBXC5), Cvoid,
        (Ptr{XCFuncType}, Float64),
        p, dens_threshold
    )
end


#
# Information
#

function Libxc_xc_func_get_info(p::Ptr{XCFuncType})
    ccall(
        (:xc_func_get_info, LIBXC5), Ptr{XCFuncInfoType},
        (Ref{XCFuncType},),
        p
    )
end


function Libxc_xc_func_info_get_family(info::Ref{XCFuncInfoType})
    ccall(
        (:xc_func_info_get_family, LIBXC5), Cint,
        (Ref{XCFuncInfoType},),
        info
    )
end
# will return an integer


function Libxc_xc_func_info_get_name(info::Ref{XCFuncInfoType})
    ccall(
        (:xc_func_info_get_name, LIBXC5), Ptr{UInt8},
        (Ref{XCFuncInfoType},),
        info
    ) |> unsafe_string
end



#
# Evaluation
#

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
    rho::AbstractVector{Float64},
    zk::AbstractVector{Float64},
    vrho::AbstractVector{Float64}
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
    rho::AbstractVector{Float64},
    sigma::AbstractVector{Float64},
    zk::AbstractVector{Float64},
    vrho::AbstractVector{Float64},
    vsigma::AbstractVector{Float64}
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
    rho::AbstractVector{Float64},
    sigma::AbstractVector{Float64},
    lapl::AbstractVector{Float64},
    tau::AbstractVector{Float64},
    vrho::AbstractVector{Float64},
    vsigma::AbstractVector{Float64},
    vlapl::AbstractVector{Float64},
    vtau::AbstractVector{Float64}
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
    rho::AbstractVector{Float64},
    sigma::AbstractVector{Float64},
    lapl::AbstractVector{Float64},
    tau::AbstractVector{Float64},
    exc::AbstractVector{Float64},
    vrho::AbstractVector{Float64},
    vsigma::AbstractVector{Float64},
    vlapl::AbstractVector{Float64},
    vtau::AbstractVector{Float64}
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