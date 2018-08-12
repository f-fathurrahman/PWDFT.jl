module PWDFT

using Printf
using LinearAlgebra
using Random
using FFTW

# constants

# CODATA: https://physics.nist.gov/cgi-bin/cuu/Value?rydhcev
const Ry2eV = 13.605693009  # Ry to eV

# CODATA: https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
# 1/bohr
const ANG2BOHR = 1.8897261254578281  # angstrom to bohr

export Ry2eV, ANG2BOHR

# Last accessed: 31st July 2018

BlochWavefunc = Array{Array{ComplexF64,2},1}
Wavefunc = Array{ComplexF64,2}
export BlochWavefunc, Wavefunc

include("Atoms.jl")
export Atoms, 
       init_atoms_xyz,
       init_atoms_xyz_string,
       get_Zatoms

include("gen_lattice_pwscf.jl")
export gen_lattice_fcc, gen_lattice_bcc, gen_lattice_trigonal, gen_lattice_trigonal_v2,
       gen_lattice_cubic, gen_lattice_bcc_v2, gen_lattice_hexagonal,
       gen_lattice_triclinic, gen_lattice_tetragonal_P, gen_lattice_tetragonal_I,
       gen_lattice_orthorhombic, gen_lattice_monoclinic,
       gen_lattice_sc, gen_lattice_rhombohedral

# FIXME Temporary workout until I found how to work with BinDeps
const SPGLIB_SO_PATH = joinpath(dirname(@__DIR__), "src", "extlibs", "spglib", "libspglib.so")
include("spglib.jl")
export spg_find_primitive,
       spg_get_ir_reciprocal_mesh
       reduce_atoms, gen_kgrid_reduced

include("KPoints.jl")
export KPoints,
       kpoints_from_file,
       kpath_from_file,
       get_special_kpoints

include("XSF_utils.jl")
export write_xsf,
       write_xsf_data3d_crystal


# General utilities
include("Utilities.jl")
export print_matrix

include("ortho_gram_schmidt.jl")
export ortho_gram_schmidt, ortho_gram_schmidt!

include("ortho_sqrt.jl")
export ortho_sqrt, ortho_sqrt!

include("ortho_check.jl")
export ortho_check

include("Ylm_real.jl")
export Ylm_real

include("Ylm_complex.jl")
export Ylm_complex

include("fft_support.jl")
export good_fft_order

#
# Plane wave basis
#
include("PWGrid.jl")
export PWGrid,
       GVectors,
       GVectorsW

#
# FFT
#
include("wrappers_fft.jl")
export R_to_G,
       G_to_R


#
# Pseudopotential
#
include("PsPot_GTH.jl")
export PsPot_GTH,
       eval_Vloc_G,
       eval_proj_G

#
# Electronic variables
#
include("Electrons.jl")
export Electrons,
       get_Zvals

include("calc_strfact.jl")
export calc_strfact

include("init_V_coulomb_G.jl")
export init_V_coulomb_G

include("calc_PspCore_ene.jl")
export calc_PspCore_ene

include("calc_E_NN.jl")
export calc_E_NN

const LIBXC_SO_PATH = joinpath(dirname(@__DIR__), "src", "extlibs", "libxc_interface.so")
include("LDA_VWN.jl")
export calc_epsxc_VWN, calc_Vxc_VWN

include("GGA_PBE.jl")
export calc_epsxc_PBE, calc_Vxc_PBE

include("Hamiltonian.jl")
export Potentials,
       Energies,
       Hamiltonian, free_electron_Hamiltonian,
       op_H, op_K, op_V_loc, op_V_Ps_loc, op_V_Ps_nloc,
       Poisson_solve,
       update!,
       PsPotNL, calc_betaNL_psi,
       calc_betaNL_psi

include("calc_rhoe.jl")
export calc_rhoe

include("Kprec.jl")
export Kprec

include("calc_energies.jl")
export calc_energies,
       calc_E_xc,
       calc_E_Hartree,
       calc_E_Ps_nloc

include("smear_FD.jl")
include("calc_Focc.jl")
include("calc_entropy.jl")
include("sum_upto_E_fermi.jl")
export smear_FD,
       calc_Focc,
       calc_entropy,
       sum_upto_E_fermi

include("calc_grad_v2.jl")
export calc_grad

#
# Diagonalization methods
#
include("diag_LOBPCG.jl")
include("diag_Emin_PCG.jl")
include("diag_davidson.jl")
export diag_LOBPCG,
       diag_Emin_PCG,
       diag_davidson

#
# Mixing functions
#
include("mix_anderson.jl")
include("mix_rpulay.jl")
export mix_anderson!,
       mix_rpulay!

include("gen_wavefunc.jl")
export rand_Wavefunc, rand_BlochWavefunc, zeros_BlochWavefunc

#
# KS solvers
#
include("KS_solve_Emin_PCG.jl")
export KS_solve_Emin_PCG!

include("KS_solve_SCF.jl")
export KS_solve_SCF!

include("KS_solve_DCM.jl")
export KS_solve_DCM!

include("KS_solve_TRDCM.jl")
export KS_solve_TRDCM!

include("CheFSI.jl")
export chebyfilt,
       get_ub_lb_lanczos

#include("do_precompile.jl")

end
