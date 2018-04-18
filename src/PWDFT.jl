__precompile__()

module PWDFT

export Ry2eV, ANG2BOHR
# constants
const Ry2eV = 13.6058         # Ry to eV
const ANG2BOHR = 1.889725989  # angstrom to bohr

export Atoms
export init_atoms_xyz
export get_Zatoms
include("Atoms.jl")

export write_xsf
export write_xsf_data3d_crystal
include("XSF_utils.jl")


#
# General utilities
#
export print_matrix
include("Utilities.jl")

export ortho_gram_schmidt
include("ortho_gram_schmidt.jl")

export ortho_check
include("ortho_check.jl")

export Ylm_real
include("Ylm_real.jl")

export Ylm_complex
include("Ylm_complex.jl")

export good_fft_order
include("fft_support.jl")

#
# Plane wave basis
#
export PWGrid
export GVectors
export GVectorsW
include("PWGrid_v03.jl")

#
# FFT
#
export c_G_to_R
export c_R_to_G
export R_to_G
export G_to_R
include("wrappers_fft.jl")

#
# Pseudopotential
#
export PsPot_GTH
export eval_Vloc_G
export eval_proj_G
include("PsPot_GTH.jl")


# Electronic variables
export Electrons
export get_Zvals
include("Electrons.jl")

export calc_strfact
include("calc_strfact.jl")

export init_V_coulomb_G
include("init_V_coulomb_G.jl")

export calc_E_NN
include("calc_E_NN.jl")

export calc_rhoe
include("calc_rhoe.jl")

const LIBXC_SO_PATH = "/home/efefer/WORKS/my_github_repos/PWDFT.jl/src/extlibs/libxc_interface.so"

export calc_epsxc_VWN, calc_Vxc_VWN
include("LDA_VWN.jl")

export calc_epsxc_PBE, calc_Vxc_PBE
include("GGA_PBE.jl")

export Potentials
export Energies
export PWHamiltonian
export op_H, op_K, op_V_loc, op_V_Ps_loc, op_V_Ps_nloc
export Poisson_solve
export update!
export PsPotNL, calc_betaNL_psi
export calc_betaNL_psi, calc_E_Ps_nloc
include("PWHamiltonian.jl")

export Kprec
include("Kprec.jl")

export calc_energies, calc_E_xc, calc_E_Hartree
include("calc_energies.jl")

export calc_grad
include("calc_grad.jl")

export diag_lobpcg
include("diag_lobpcg.jl")

export diag_Emin_PCG
include("diag_Emin_PCG.jl")

export andersonmix!
include("andersonmix.jl")

export KS_solve_Emin_PCG!
include("KS_solve_Emin_PCG.jl")

export KS_solve_SCF!
include("KS_solve_SCF.jl")

export chebyfilt
export get_ub_lb_lanczos
include("CheFSI.jl")

export KS_solve_DCM!
include("KS_solve_DCM.jl")


end
