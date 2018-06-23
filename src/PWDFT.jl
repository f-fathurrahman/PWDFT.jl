#
# Going from 0.6 to 0.7
#
# Void -> Nothing
# Complex128 -> ComplexF64
# find! -> findall!
# FFT Plan types

__precompile__()

module PWDFT

using Printf
using LinearAlgebra
using Random
using FFTW

export Ry2eV, ANG2BOHR
# constants
const Ry2eV = 13.6058         # Ry to eV
const ANG2BOHR = 1.889725989  # angstrom to bohr

export Atoms
export init_atoms_xyz
export init_atoms_xyz_string
export get_Zatoms
include("Atoms.jl")

export gen_lattice_fcc, gen_lattice_bcc, gen_lattice_trigonal, gen_lattice_trigonal_v2,
       gen_lattice_cubic, gen_lattice_bcc_v2, gen_lattice_hexagonal,
       gen_lattice_triclinic, gen_lattice_tetragonal_P, gen_lattice_tetragonal_I,
       gen_lattice_orthorhombic, gen_lattice_monoclinic,
       gen_lattice_sc
include("gen_lattice_pwscf.jl")

# FIXME Temporary workout until I found how to work with BinDeps
const SPGLIB_SO_PATH = "/home/efefer/WORKS/my_github_repos/PWDFT.jl/src/extlibs/spglib/libspglib.so"
export spg_find_primitive
export spg_get_ir_reciprocal_mesh
export reduce_atoms, gen_kgrid_reduced
include("spglib.jl")


export KPoints
export kpoints_from_file
export kpath_from_file
export get_special_kpoints
include("KPoints.jl")


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
include("PWGrid.jl")

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

#
# Electronic variables
#
export Electrons
export get_Zvals
include("Electrons.jl")

export calc_strfact
include("calc_strfact.jl")

export init_V_coulomb_G
include("init_V_coulomb_G.jl")

export calc_PspCore_ene
include("calc_PspCore_ene.jl")

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
export calc_betaNL_psi
include("PWHamiltonian.jl")

export Kprec
include("Kprec.jl")

export calc_energies, calc_E_xc, calc_E_Hartree, calc_E_Ps_nloc
include("calc_energies.jl")

export smear_FD, calc_Focc, calc_entropy, sum_upto_E_fermi
include("smear_FD.jl")
include("calc_Focc.jl")
include("calc_entropy.jl")
include("sum_upto_E_fermi.jl")

export calc_grad
include("calc_grad_v2.jl")

#
# Diagonalization methods
#
export diag_lobpcg
include("diag_lobpcg.jl")

export diag_Emin_PCG
include("diag_Emin_PCG.jl")

export diag_davidson
include("diag_davidson.jl")

#
# Mixing functions
#
export andersonmix!
include("andersonmix.jl")

#
# KS solvers
#
export KS_solve_Emin_PCG!
include("KS_solve_Emin_PCG.jl")

export KS_solve_SCF!
include("KS_solve_SCF.jl")

export KS_solve_SCF_smearing!
include("KS_solve_SCF_smearing.jl")

export KS_solve_DCM!
include("KS_solve_DCM.jl")

export chebyfilt
export get_ub_lb_lanczos
include("CheFSI.jl")

#include("do_precompile.jl")

end
