module PWDFT2

using Printf
using LinearAlgebra
using Random
using FFTW

include("../extlibs/extlibs.jl") # Load library dependencies

# constants

# CODATA: https://physics.nist.gov/cgi-bin/cuu/Value?rydhcev
const Ry2eV = 13.605693009  # Ry to eV

# CODATA: https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
# 1/bohr
const ANG2BOHR = 1.8897261254578281  # angstrom to bohr

export Ry2eV, ANG2BOHR

# Last accessed: 31st July 2018

const BlochWavefunc = Array{Array{ComplexF64,2},1}
const Wavefunc = Array{ComplexF64,2}
export BlochWavefunc, Wavefunc

include("../src/Atoms.jl")
export Atoms, 
       init_atoms_xyz,
       init_atoms_xyz_string,
       get_Zatoms

include("../src/gen_lattice_pwscf.jl")
export gen_lattice_fcc, gen_lattice_bcc, gen_lattice_trigonal, gen_lattice_trigonal_v2,
       gen_lattice_cubic, gen_lattice_bcc_v2, gen_lattice_hexagonal,
       gen_lattice_triclinic, gen_lattice_tetragonal_P, gen_lattice_tetragonal_I,
       gen_lattice_orthorhombic, gen_lattice_monoclinic,
       gen_lattice_sc, gen_lattice_rhombohedral

include("../src/spglib.jl")
export spg_find_primitive,
       spg_get_ir_reciprocal_mesh,
       reduce_atoms, gen_kgrid_reduced

include("../src/KPoints.jl")
export KPoints,
       kpoints_from_file,
       kpath_from_file,
       get_special_kpoints

include("../src/XSF_utils.jl")
export write_xsf,
       write_xsf_data3d_crystal


# General utilities
include("../src/Utilities.jl")
export print_matrix

include("../src/ortho_gram_schmidt.jl")
export ortho_gram_schmidt, ortho_gram_schmidt!

include("../src/ortho_sqrt.jl")
export ortho_sqrt, ortho_sqrt!

include("../src/ortho_check.jl")
export ortho_check

include("../src/Ylm_real.jl")
export Ylm_real

include("../src/Ylm_complex.jl")
export Ylm_complex

include("../src/fft_support.jl")
export good_fft_order

#
# Plane wave basis
#
include("../src/PWGrid.jl")
export PWGrid,
       GVectors,
       GVectorsW,
       op_nabla, op_nabla_dot

#
# FFT
#
include("../src/wrappers_fft.jl")
export R_to_G,
       G_to_R


include("PsPot_GTH_mod.jl")
export PsPot_GTH,
       eval_Vloc_G,
       eval_proj_G

#
# Electronic variables
#
include("../src/Electrons.jl")
export Electrons,
       get_Zvals

include("../src/calc_strfact.jl")
export calc_strfact

include("../src/init_V_coulomb_G.jl")
export init_V_coulomb_G

include("../src/calc_PspCore_ene.jl")
export calc_PspCore_ene

include("../src/calc_E_NN.jl")
export calc_E_NN

mutable struct XCFuncType
end

include("../src/LDA_VWN.jl")
export calc_epsxc_VWN, calc_Vxc_VWN

include("../src/GGA_PBE.jl")
export calc_epsxc_PBE, calc_Vxc_PBE

include("Energies_mod.jl")
export Energies

include("Potentials_mod.jl")
export Energies

include("../src/PsPotNL.jl")
export PsPotNL, calc_betaNL_psi

include("../src/free_electron_Hamiltonian.jl")
export free_electron_Hamiltonian

include("../src/Hamiltonian.jl")
export Hamiltonian, update!

include("../src/op_K.jl")
include("../src/op_V_loc.jl")
include("../src/op_V_Ps_nloc.jl")
include("../src/op_H.jl")
export op_H, op_K, op_V_loc, op_V_Ps_loc, op_V_Ps_nloc

include("../src/Poisson_solve.jl")
export Poisson_solve

include("../src/calc_rhoe.jl")
export calc_rhoe

include("../src/Kprec.jl")
export Kprec

include("../src/precKerker.jl")
export precKerker

include("../src/calc_energies.jl")
export calc_energies,
       calc_E_xc,
       calc_E_Hartree,
       calc_E_Ps_nloc

include("../src/smear_FD.jl")
include("../src/calc_Focc.jl")
include("../src/calc_entropy.jl")
include("../src/sum_upto_E_fermi.jl")
export smear_FD,
       calc_Focc,
       calc_entropy,
       sum_upto_E_fermi

include("../src/calc_grad_v2.jl")
export calc_grad

#
# Diagonalization methods
#
include("../src/diag_LOBPCG.jl")
include("../src/diag_Emin_PCG.jl")
include("../src/diag_davidson.jl")
export diag_LOBPCG, diag_LOBPCG!,
       diag_Emin_PCG, diag_Emin_PCG!,
       diag_davidson, diag_davidson!

#
# Mixing functions
#
include("../src/mix_anderson.jl")
include("../src/mix_rpulay.jl")
export mix_anderson!,
       mix_rpulay!

include("../src/gen_wavefunc.jl")
export rand_Wavefunc, rand_BlochWavefunc, zeros_BlochWavefunc

#
# KS solvers
#
include("../src/KS_solve_Emin_PCG.jl")
export KS_solve_Emin_PCG!

include("../src/KS_solve_SCF.jl")
export KS_solve_SCF!

include("../src/KS_solve_DCM.jl")
export KS_solve_DCM!

include("../src/KS_solve_TRDCM.jl")
export KS_solve_TRDCM!

include("../src/CheFSI.jl")
export chebyfilt,
       get_ub_lb_lanczos

#include("do_precompile.jl")

end
