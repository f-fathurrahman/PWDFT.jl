module PWDFT_cuda

using Printf
using FFTW
using LinearAlgebra

using CUDAnative
using CuArrays

using PWDFT

import CuArrays: allowscalar, @allowscalar
allowscalar(false)

include("utils_CuArrays.jl")

include("kernels_fft_grid.jl")

include("CuElectrons.jl")
export CuElectrons

include("CuPWGrid.jl")
export CuPWGrid, CuGVectors, CuGVectorsW

include("cu_calc_strfact.jl")

include("cu_wrappers_fft.jl")

include("cu_types_aliases.jl")
export CuWavefunc, CuBlochWavefunc

include("cu_ortho_gram_schmidt.jl")
include("cu_ortho_check.jl")

include("cu_XC_functionals_internal.jl")
include("cu_LDA_VWN_internal.jl")
include("cu_GGA_PBE_internal.jl")

include("CuPotentials.jl")
export CuPotentials

include("CuPsPotNL.jl")
export CuPsPotNL

include("CuHamiltonian.jl")
export CuHamiltonian

include("cu_gen_wavefunc.jl")
export zeros_CuBlochWavefunc, rand_CuBlochWavefunc

include("cu_Kprec.jl")
include("cu_op_K.jl")
include("cu_op_V_loc.jl")
include("cu_op_V_Ps_nloc.jl")
include("cu_op_H.jl")

include("cu_calc_rhoe.jl")

include("cu_Poisson_solve.jl")

include("cu_calc_energies.jl")

include("cu_calc_grad.jl")

include("cu_KS_solve_Emin_PCG.jl")

end  # end of module

