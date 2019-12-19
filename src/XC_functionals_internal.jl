# Adapted from PWSCF
# I found slight difference with Libxc for PBE correlation.

include("XC_funcs/XC_c_pbe.jl")
include("XC_funcs/XC_c_pbe_E.jl")

include("XC_funcs/XC_c_pbe_spin.jl")
include("XC_funcs/XC_c_pbe_spin_E.jl")

include("XC_funcs/XC_c_pw.jl")
include("XC_funcs/XC_c_pw_E.jl")

include("XC_funcs/XC_c_pw_spin.jl")
include("XC_funcs/XC_c_pw_spin_E.jl")

include("XC_funcs/XC_c_vwn.jl")

include("XC_funcs/XC_c_vwn_spin.jl")


include("XC_funcs/XC_x_pbe.jl")
include("XC_funcs/XC_x_pbe_E.jl")

include("XC_funcs/XC_x_slater.jl")

include("XC_funcs/XC_x_slater_spin.jl")
include("XC_funcs/XC_x_slater_spin_E.jl")
