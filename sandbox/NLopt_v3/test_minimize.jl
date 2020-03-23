using LinearAlgebra
using Printf
using Random

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("calc_energies_grad.jl")
include("create_Ham.jl")
include("KS_solve_Emin_PCG_new.jl")
include("linmin_quad.jl")
include("linmin_grad.jl")
include("linmin_debug.jl")

function main()
    Random.seed!(1234)
    
    #Ham = create_Ham_H2()
    #Ham = create_Ham_H_atom()
    Ham = create_Ham_Si_fcc()
    #Ham = create_Ham_GaAs()
    #Ham = create_Ham_NH3()
    #Ham = create_Ham_ZnO()

    psiks = rand_BlochWavefunc( Ham )
    KS_solve_Emin_PCG_new!( Ham, psiks, startingrhoe=:gaussian, skip_initial_diag=false, NiterMax=100 )
    #KS_solve_Emin_PCG_new!( Ham, psiks, startingrhoe=:random, skip_initial_diag=true, NiterMax=100 )
    #linmin_debug!( Ham, psiks, startingrhoe=:random, skip_initial_diag=true )
    
    #KS_solve_SCF!( Ham, mix_method="anderson" )
    #KS_solve_Emin_PCG!( Ham )
end

main()
