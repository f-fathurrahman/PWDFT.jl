using LinearAlgebra
using Printf
using Random

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include(joinpath(DIR_PWDFT, "sandbox", "KS_solve_SCF_NLsolve.jl"))
include(joinpath(DIR_PWDFT, "utilities", "PWSCF.jl"))

include("MinimizeParams.jl")
include("calc_energies_grad.jl")
include("create_Ham.jl")
include("KS_solve_Emin_PCG_new.jl")
include("linmin_quad.jl")
include("linmin_grad.jl")
include("linmin_debug.jl")
include("linmin_armijo.jl")
include("KS_solve_Emin_PCG_dot.jl")

function main()
    Random.seed!(1234)
    #Random.seed!(1111)

    #Ham = create_Ham_H2()
    #Ham = create_Ham_H_atom()
    Ham = create_Ham_Si_fcc()
    #Ham = create_Ham_GaAs()
    #Ham = create_Ham_NH3()
    #Ham = create_Ham_ZnO()

    #println(Ham)

    psiks = rand_BlochWavefunc( Ham )
    #KS_solve_Emin_PCG_new!( Ham, psiks, NiterMax=100 ) #, etot_conv_thr=1e-7 )
    #KS_solve_Emin_PCG_new!( Ham, psiks, startingrhoe=:random, skip_initial_diag=true, NiterMax=150 )
    #linmin_debug!( Ham, psiks, startingrhoe=:random, skip_initial_diag=true )
    
    KS_solve_SCF!( Ham, mix_method="anderson", betamix=0.5 , etot_conv_thr=1e-8 )
    #KS_solve_SCF_NLsolve!( Ham )
    #KS_solve_Emin_PCG!( Ham, startingrhoe=:random, skip_initial_diag=true, i_cg_beta=4 )
    #KS_solve_Emin_PCG!( Ham, i_cg_beta=2 )
    #KS_solve_Emin_PCG_dot!( Ham, psiks, startingrhoe=:random, skip_initial_diag=true, etot_conv_thr=1e-8 )
    #KS_solve_Emin_PCG_dot!( Ham, psiks, etot_conv_thr=1e-5 )
end

main()