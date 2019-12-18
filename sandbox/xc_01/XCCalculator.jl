abstract type AbstractXCCalculator
end

struct XCCalculator <: AbstractXCCalculator
end

struct LibxcXCCalculator <: AbstractXCCalculator
end

include("XC_x_slater.jl")
include("XC_c_pw.jl")
include("XC_c_vwn.jl")

include("XC_x_slater_spin.jl")
include("XC_c_pw_spin.jl")
include("XC_c_vwn_spin.jl")

include("XC_x_pbe.jl")
include("XC_c_pbe.jl")

include("XC_c_pbe_spin.jl")