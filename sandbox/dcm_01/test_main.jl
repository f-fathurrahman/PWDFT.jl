using Printf
using LinearAlgebra
using Random
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("create_Ham.jl")

include("KS_solve_DCM_01.jl")
include("KS_solve_TRDCM_01.jl")
include("my_dcm.jl")

function main()

    Random.seed!(1234)

    #Ham = create_Ham_Si_fcc()
    Ham = create_Ham_CO()
    #Ham = create_Ham_H2()
    #Ham = create_Ham_GaAs_v1()

    #my_dcm!( Ham )

    #KS_solve_DCM_01!( Ham, NiterMax=30, MaxInnerSCF=3 )
    #@time KS_solve_DCM_01!( Ham, NiterMax=20, MaxInnerSCF=3 )

    #KS_solve_TRDCM_01!( Ham, NiterMax=20 )
    #KS_solve_Emin_PCG!( Ham )
    KS_solve_SCF!( Ham, mix_method="anderson" )

end

main()