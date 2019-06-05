using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("mix_broyden.jl")
include("create_Ham.jl")
include("KS_solve_SCF_broyden.jl")

function test_main()

    #Ham = create_Ham_H2()
    Ham = create_Ham_Si_fcc()
    #Ham = create_Ham_Si_fcc(xcfunc="PBE")
    #Ham = create_Ham_GaAs_v1()
    #Ham = create_Ham_GaAs_v2()
    #Ham = create_Ham_CO()
    println(Ham)

    Random.seed!(1234)
    
    KS_solve_SCF_broyden!(Ham, mix_method="rpulay", mixdim=8, betamix=0.7)
    #KS_solve_Emin_PCG!(Ham)

end


test_main()
