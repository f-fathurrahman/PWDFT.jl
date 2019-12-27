# Adapted from PWSCF
# I found slight difference with Libxc for PBE correlation.

include("cu_XC_funcs/cu_XC_c_pbe.jl")
include("cu_XC_funcs/cu_XC_c_pbe_E.jl")

include("cu_XC_funcs/cu_XC_c_pbe_spin.jl")
include("cu_XC_funcs/cu_XC_c_pbe_spin_E.jl")

include("cu_XC_funcs/cu_XC_c_pw.jl")
include("cu_XC_funcs/cu_XC_c_pw_E.jl")

include("cu_XC_funcs/cu_XC_c_pw_spin.jl")
include("cu_XC_funcs/cu_XC_c_pw_spin_E.jl")

include("cu_XC_funcs/cu_XC_c_vwn.jl")

include("cu_XC_funcs/cu_XC_c_vwn_spin.jl")

include("cu_XC_funcs/cu_XC_x_pbe.jl")
include("cu_XC_funcs/cu_XC_x_pbe_E.jl")

include("cu_XC_funcs/cu_XC_x_slater.jl")

include("cu_XC_funcs/cu_XC_x_slater_spin.jl")
include("cu_XC_funcs/cu_XC_x_slater_spin_E.jl")

