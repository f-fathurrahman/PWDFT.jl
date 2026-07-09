struct NoneXCCalculator <: AbstractXCCalculator
    x_id::Int64
    c_id::Int64
    family::Symbol
    need_gradient::Bool
end

function NoneXCCalculator()
    return NoneXCCalculator(
        0, 0, :none, false
    )
end

# In-place version, both epsxc and Vxc, spinpol
function calc_epsxc_Vxc_VWN!(
    xc_calc::NoneXCCalculator,
    Rhoe::Array{Float64,2},
    epsxc::Vector{Float64},
    Vxc::Array{Float64,2}
)
    fill!(epsxc, 0.0)
    fill!(Vxc, 0.0)
    return
end


function calc_epsxc_Vxc_VWN!(
    xc_calc::NoneXCCalculator,
    Rhoe::AbstractVector{Float64},
    epsxc::AbstractVector{Float64},
    Vxc::AbstractVector{Float64}
)
    fill!(epsxc, 0.0)
    fill!(Vxc, 0.0)
    return
end


function calc_epsxc_Vxc_VWN(
    xc_calc::NoneXCCalculator,
    Rhoe::AbstractVector{Float64}
)
    Npoints = size(Rhoe, 1)
    epsxc = zeros(Float64, Npoints)
    Vxc = zeros(Float64, Npoints)
    return epsxc, Vxc
end

#
# epsxc only
#
function calc_epsxc_VWN!(
    xc_calc::NoneXCCalculator,
    Rhoe::AbstractVector{Float64},
    epsxc::AbstractVector{Float64}
)
    fill!(epsxc, 0.0)
    return
end

function calc_epsxc_VWN(
    xc_calc::NoneXCCalculator,
    Rhoe
)
    epsxc = zeros(Float64, size(Rhoe))
    return epsxc
end


function calc_epsxc_VWN!(
    xc_calc::NoneXCCalculator,
    Rhoe::Matrix{Float64},
    epsxc::Vector{Float64}
)
    fill!(epsxc, 0.0)
    return
end


function calc_Vxc_VWN!(
    xc_calc::NoneXCCalculator,
    Rhoe::AbstractVector{Float64},
    Vxc::AbstractVector{Float64}
)
    fill!(Vxc, 0.0)
    return
end


function calc_Vxc_VWN(
    xc_calc::NoneXCCalculator,
    Rhoe::AbstractVector{Float64}
)
    Vxc = zeros(Float64, size(Rhoe))
    return Vxc
end


function calc_Vxc_VWN!(
    xc_calc::NoneXCCalculator,
    Rhoe::Array{Float64,2},
    Vxc::Array{Float64,2}
)
    fill!(Vxc, 0.0)
    return
end
