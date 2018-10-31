module PWDFT

using Printf
using LinearAlgebra
using Random
using FFTW

# Load library dependencies
# NOTE: This file must be edited manually
include("../extlibs/extlibs.jl") 

# constants
#
# CODATA: https://physics.nist.gov/cgi-bin/cuu/Value?rydhcev
const Ry2eV = 13.605693009  # Ry to eV
#
# CODATA: https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
# 1/bohr
const ANG2BOHR = 1.8897261254578281  # angstrom to bohr

export Ry2eV, ANG2BOHR

# Some aliases
const BlochWavefunc = Array{Array{ComplexF64,2},1}
const Wavefunc = Array{ComplexF64,2}
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

include("spglib.jl")
export spg_find_primitive,
       spg_get_ir_reciprocal_mesh,
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
       GVectorsW,
       op_nabla, op_nabla_dot

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
       print_ebands,
       get_Zvals

include("calc_strfact.jl")
export calc_strfact

include("init_V_coulomb_G.jl")
export init_V_coulomb_G

include("calc_PspCore_ene.jl")
export calc_PspCore_ene

include("calc_E_NN.jl")
export calc_E_NN

mutable struct XCFuncType
end

export XCFuncType

include("LDA_VWN.jl")
export calc_epsxc_VWN, calc_Vxc_VWN

include("GGA_PBE.jl")
export calc_epsxc_PBE, calc_Vxc_PBE

include("Energies.jl")
export Energies

include("Potentials.jl")
export Energies

include("PsPotNL.jl")
export PsPotNL, calc_betaNL_psi

include("free_electron_Hamiltonian.jl")
export free_electron_Hamiltonian

include("Hamiltonian.jl")
export Hamiltonian, update!

include("op_K.jl")
include("op_V_loc.jl")
include("op_V_Ps_nloc.jl")
include("op_H.jl")
export op_H, op_K, op_V_loc, op_V_Ps_loc, op_V_Ps_nloc

include("Poisson_solve.jl")
export Poisson_solve

include("calc_rhoe.jl")
export calc_rhoe

include("Kprec.jl")
export Kprec

include("precKerker.jl")
export precKerker

include("calc_energies.jl")
export calc_energies,
       calc_E_xc,
       calc_E_Hartree,
       calc_E_Ps_nloc

include("calc_Focc.jl")
include("occupations.jl")
include("calc_entropy.jl")
include("sum_upto_E_fermi.jl")
export calc_Focc,
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
export diag_LOBPCG, diag_LOBPCG!,
       diag_Emin_PCG, diag_Emin_PCG!,
       diag_davidson, diag_davidson!

#
# Mixing functions
#
include("mix_anderson.jl")
include("mix_pulay.jl")
include("mix_rpulay.jl")
include("mix_ppulay.jl")
export mix_anderson!,
       mix_pulay!,
       mix_rpulay!,
       mix_ppulay!

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

end
