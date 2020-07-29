using Printf
using LinearAlgebra
using Random
using FFTW

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("PWGridGamma.jl")
include("wrappers_fft_gamma.jl")
include("ortho_GS_gamma.jl")
include("PsPotNLGamma.jl")
include("HamiltonianGamma.jl")
include("BlochWavefuncGamma.jl")
include("calc_rhoe_gamma.jl")
include("Poisson_solve_gamma.jl")
include("op_K_gamma.jl")
include("op_V_loc_gamma.jl")
include("op_V_Ps_nloc_gamma.jl")
include("op_H_gamma.jl")
include("calc_energies_gamma.jl")
include("calc_grad_gamma.jl")

include("setup_guess_wavefunc.jl")

include("KS_solve_Emin_PCG_dot.jl")
include("calc_energies_grad.jl")
include("linmin_grad.jl")

include("KS_solve_Emin_PCG_dot_gamma.jl")
include("calc_energies_grad_gamma.jl")
include("linmin_grad_gamma.jl")

include("unfold_BlochWavefuncGamma.jl")

include(joinpath(DIR_PWDFT, "sandbox", "get_default_psp.jl"))

function main(filename)

    Random.seed!(1234)

    atoms = Atoms(ext_xyz_file=filename)
    pspfiles = get_default_psp(atoms)
    
    ecutwfc = 15.0
    
    Ham = HamiltonianGamma(atoms, pspfiles, ecutwfc )
    psis = randn_BlochWavefuncGamma(Ham)
    @time KS_solve_Emin_PCG_dot!( Ham, psis, NiterMax=200 )
    
    #Ham_ = Hamiltonian( atoms, pspfiles, ecutwfc, use_symmetry=false )
    #psiks = unfold_BlochWavefuncGamma( Ham.pw, Ham_.pw, psis )
    #@time KS_solve_Emin_PCG_dot!( Ham_, psiks, startingrhoe=:random, skip_initial_diag=true, NiterMax=200 )
    #@time KS_solve_Emin_PCG!( Ham_, psiks, startingrhoe=:random, skip_initial_diag=true )
    #@time KS_solve_Emin_PCG!( Ham_, psiks )

end

#main(joinpath(DIR_STRUCTURES, "DATA_G2_mols", "H2O.xyz"))
#main(joinpath(DIR_STRUCTURES, "DATA_G2_mols", "NH3.xyz"))
#main(joinpath(DIR_STRUCTURES, "DATA_G2_mols", "SiH4.xyz"))
#main(joinpath(DIR_STRUCTURES, "DATA_G2_mols", "N2H4.xyz"))
main(joinpath(DIR_STRUCTURES, "DATA_G2_mols", "CO2.xyz"))