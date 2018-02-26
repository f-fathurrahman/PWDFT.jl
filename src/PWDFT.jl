__precompile__()

module PWDFT

export Ry2eV, ANG2BOHR
# constants
const Ry2eV = 13.6058         # Ry to eV
const ANG2BOHR = 1.889725989  # angstrom to bohr

export Atoms
export init_atoms_xyz
include("Atoms.jl")


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


#
# Plane wave basis
#
export PWGrid
export GVectors
export GVectorsW
include("PWGrid_v03.jl")

#
# FFTW
#
export c_G_to_R
export c_R_to_G
export R_to_G
export G_to_R
include("wrappers_fft.jl")

#
# Pseudopotential
#
export PsPot_HGH
export PsPot_GTH
export eval_Vloc_G
export eval_proj_G

include("PsPot_HGH.jl")
include("PsPot_GTH.jl")


export calc_strfact
include("calc_strfact.jl")

export init_V_coulomb_G
include("init_V_coulomb_G.jl")

export calc_E_NN
include("calc_E_NN.jl")

export calc_rhoe
include("calc_rhoe.jl")

export PotentialsT
export EnergiesT
export PWHamiltonian
export op_H, op_K, op_V_loc, op_V_Ps_loc
export Poisson_solve
export update!
include("PWHamiltonian.jl")

export Kprec
include("Kprec.jl")

export calc_energies
include("calc_energies.jl")

export calc_grad
include("calc_grad.jl")

export diag_lobpcg
include("diag_lobpcg.jl")

export diag_Emin_PCG
include("diag_Emin_PCG.jl")

export KS_solve_Emin_PCG!
include("KS_solve_Emin_PCG.jl")

export KS_solve_SCF!
include("KS_solve_SCF.jl")

end
