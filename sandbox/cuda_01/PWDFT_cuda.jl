using Printf
using FFTW
using LinearAlgebra

using CUDAnative
using CuArrays

using PWDFT

import CuArrays: allowscalar
allowscalar(false)

include("CuPWGrid.jl")

include("cu_calc_strfact.jl")

include("cu_wrappers_fft.jl")

include("cu_types_aliases.jl")

include("cu_ortho_gram_schmidt.jl")
include("cu_gen_wavefunc.jl")
include("cu_ortho_check.jl")

include("CuPotentials.jl")
include("CuHamiltonian.jl")
