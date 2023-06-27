module PWDFT

using Printf
using LinearAlgebra
using Random
using FFTW
using LightXML
using OffsetArrays

# constants
#
# CODATA: https://physics.nist.gov/cgi-bin/cuu/Value?rydhcev
const Ry2eV = 13.605693009  # Ry to eV
const Ha2eV = 2*Ry2eV
const eV2Ha = 1/Ha2eV

#
# CODATA: https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
# 1/bohr
const ANG2BOHR = 1.8897261254578281  # angstrom to bohr
const BOHR2ANG = 1/ANG2BOHR

export Ry2eV, Ha2eV, eV2Ha, ANG2BOHR, BOHR2ANG

# Some aliases
const BlochWavefunc = Array{Array{ComplexF64,2},1}
const Wavefunc = Array{ComplexF64,2}
export BlochWavefunc, Wavefunc

include("Atoms.jl")
export Atoms,
       init_atoms_xyz,
       init_atoms_xyz_string,
       get_Zatoms

include("gen_lattice.jl")
export gen_lattice_fcc, gen_lattice_bcc, gen_lattice_trigonal, gen_lattice_trigonal_v2,
       gen_lattice_cubic, gen_lattice_bcc_v2, gen_lattice_hexagonal,
       gen_lattice_triclinic, gen_lattice_tetragonal_P, gen_lattice_tetragonal_I,
       gen_lattice_orthorhombic, gen_lattice_monoclinic,
       gen_lattice_sc, gen_lattice_rhombohedral

include("SymmetryInfo.jl")
export SymmetryInfo

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

include("Ylm_real.jl")
export Ylm_real

include("Ylm_complex.jl")
export Ylm_complex

include("Ylm_real_qe.jl")
export Ylm_real_qe!

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

include("gamma_only/PWGridGamma.jl")
export PWGridGamma,
       GVectorsGamma,
       GVectorsWGamma

#
# FFT
#
include("wrappers_fft.jl")
export R_to_G, R_to_G!,
       G_to_R, G_to_R!

#
# Pseudopotential
#
abstract type AbstractPsPot end

include("PsPot_GTH.jl")
export PsPot_GTH,
       eval_Vloc_G, eval_Vloc_G!,
       eval_proj_G,
       write_psp10

include("PAWData_UPF.jl")
export PAWData_UPF

include("PsPot_UPF.jl")
export PsPot_UPF, calc_Natomwfc

export is_using_extension_gth, is_using_extension_upf

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

include("XCCalculator.jl")
export AbstractXCCalculator, XCCalculator, LibxcXCCalculator

include("Libxc_old.jl")

include("LDA_VWN.jl")
export calc_epsxc_VWN, calc_Vxc_VWN
export calc_epsxc_Vxc_VWN, calc_epsxc_Vxc_VWN!

include("GGA_PBE.jl")
export calc_epsxc_PBE, calc_Vxc_PBE

include("XC_functionals_internal.jl")
include("LDA_VWN_internal.jl")
include("GGA_PBE_internal.jl")

include("Energies.jl")
export Energies

include("Potentials.jl")
export Potentials

include("PsPotNL.jl")
export PsPotNL, calc_betaNL_psi

include("PAWAtomicSphere.jl")
include("PAWVariables.jl")
export PAWAtomicSphere, PAWVariables

include("qvan2.jl")
include("PsPotNL_UPF.jl")
export PsPotNL_UPF, qvan2!

include("gamma_only/PsPotNLGamma.jl")
export PsPotNLGamma

include("RhoeSymmetrizer.jl")
export RhoeSymmetrizer,
       symmetrize_rhoe!

include("calc_rhoe_core.jl")
export calc_rhoe_core!

include("Hamiltonian.jl")
export Hamiltonian, HamiltonianOptions, update!

include("gamma_only/HamiltonianGamma.jl")
export HamiltonianGamma

include("gamma_only/BlochWavefuncGamma.jl")
export BlochWavefuncGamma,
       rand_BlochWavefuncGamma,
       randn_BlochWavefuncGamma,
       dot_gamma, overlap_gamma,
       ortho_sqrt_gamma, ortho_sqrt_gamma!,
       ortho_check_gamma

include("gamma_only/unfold_BlochWavefuncGamma.jl")
export unfold_BlochWavefuncGamma

include("ortho_sqrt.jl")
export ortho_sqrt, ortho_sqrt!

include("ortho_check.jl")
export ortho_check

# We need Hamiltonian to be declared here (?)
include("MGGA_SCAN.jl")
export calc_epsxc_SCAN, calc_Vxc_SCAN!
export calc_KEdens!

include("op_K.jl")
include("op_V_loc.jl")
include("op_V_Ps_nloc.jl")
include("op_Vtau.jl")
include("op_H.jl")
export op_H, op_K, op_V_loc, op_V_Ps_loc, op_V_Ps_nloc,
       op_H!, op_K!, op_V_loc!, op_V_Ps_loc!, op_V_Ps_nloc!,
       op_Vtau, op_Vtau!

include("gamma_only/op_K_gamma.jl")
include("gamma_only/op_V_loc_gamma.jl")
include("gamma_only/op_V_Ps_nloc_gamma.jl")
include("gamma_only/op_H_gamma.jl")

include("PAW_symmetrize.jl")
export PAW_symmetrize!

include("PAW_symmetrize_ddd.jl")
export PAW_symmetrize_ddd!

#
# Stuff related to PAW_potential
#
include("PAW_atomic_becsum.jl")
export PAW_atomic_becsum, PAW_atomic_becsum!

include("PAW_rho_lm.jl")
export PAW_rho_lm!

include("radial_hartree.jl")
include("PAW_h_potential.jl")
export PAW_h_potential!

include("PAW_lm2rad.jl")
include("PAW_rad2lm.jl")
export PAW_lm2rad!, PAW_rad2lm!

include("PAW_xc_potential.jl")
export PAW_xc_potential!

include("PAW_potential.jl")
export PAW_potential!


include("Poisson_solve.jl")
export Poisson_solve, Poisson_solve!

include("gamma_only/Poisson_solve_gamma.jl")


include("calc_rhoe.jl")
export calc_rhoe, calc_rhoe!

include("gamma_only/calc_rhoe_gamma.jl")

include("Kprec.jl")
export Kprec, Kprec!, Kprec_inplace!

include("calc_energies.jl")
export calc_energies,
       calc_E_kin,
       calc_E_local,
       calc_E_Ps_nloc
# calc_E_xc, calc_E_Hartree # Not needed?

include("gamma_only/calc_energies_gamma.jl")


include("occupations.jl")
export calc_Focc,
       calc_entropy,
       sum_upto_E_fermi

include("calc_grad.jl")
export calc_grad

include("gamma_only/calc_grad_gamma.jl")
export calc_grad! # FIXME, also implement for kpt-version

#
# Diagonalization methods
#
include("diag_LOBPCG.jl")
include("diag_Emin_PCG.jl")
include("diag_davidson.jl")
export diag_LOBPCG, diag_LOBPCG!,
       diag_Emin_PCG, diag_Emin_PCG!,
       diag_davidson, diag_davidson!

include("gamma_only/diag_Emin_PCG_gamma.jl")


#
# Mixing functions
#
include("mix_anderson.jl")
include("mix_adaptive.jl")
include("mix_pulay.jl")
include("mix_rpulay.jl")
include("mix_ppulay.jl")
include("mix_broyden.jl")
export mix_anderson!,
       mix_adaptive!,
       mix_pulay!,
       mix_rpulay!,
       mix_ppulay!,
       mix_broyden!

# Newer
include("BroydenMixer.jl")
export do_mix!,
       BroydenMixer



include("gen_wavefunc.jl")
export rand_Wavefunc, rand_BlochWavefunc, zeros_BlochWavefunc

include("guess_rhoe.jl")
export guess_rhoe, guess_rhoe_atomic

#
# KS solvers
#
include("KS_solve_Emin_PCG.jl")
export KS_solve_Emin_PCG!

include("gamma_only/setup_guess_wavefunc.jl")
include("gamma_only/linmin_grad_gamma.jl")
include("gamma_only/calc_energies_grad_gamma.jl")
include("gamma_only/KS_solve_Emin_PCG_dot_gamma.jl")
export KS_solve_Emin_PCG_dot!,
       setup_guess_wavefunc!

include("KS_solve_SCF.jl")
export KS_solve_SCF!

include("gamma_only/KS_solve_SCF_gamma.jl")

include("KS_solve_SCF_potmix.jl")
export KS_solve_SCF_potmix!

include("gamma_only/KS_solve_SCF_potmix_gamma.jl")

include("KS_solve_DCM.jl")
export KS_solve_DCM!

include("KS_solve_TRDCM.jl")
export KS_solve_TRDCM!

include("CheFSI.jl")
export chebyfilt,
       get_ub_lb_lanczos

include("read_psiks.jl")
export read_psiks


include("electrons_scf.jl")
export electrons_scf!
export atomic_rho_g
export update_from_rhoe!
export op_S, op_S!


include("calc_forces_NN.jl")
include("calc_forces_Ps_loc.jl")
include("calc_forces_Ps_nloc.jl")
include("calc_forces.jl")
export calc_forces_NN, calc_forces_NN!
export calc_forces_Ps_loc, calc_forces_Ps_loc!
export calc_forces_Ps_nloc, calc_forces_Ps_nloc!
export calc_forces

include("gamma_only/calc_forces_NN_gamma.jl")
include("gamma_only/calc_forces_Ps_loc_gamma.jl")
include("gamma_only/calc_forces_Ps_nloc_gamma.jl")
include("gamma_only/calc_forces_gamma.jl")


# Some utilities
include("get_default_psp.jl")
export get_default_PsPot_GTH

end
