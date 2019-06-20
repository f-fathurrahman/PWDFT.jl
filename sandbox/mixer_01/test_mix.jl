using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("precKerker.jl")
include("create_Ham.jl")
include("mix_adaptive.jl")
include("KS_solve_SCF_rhomix_v2.jl")
include("KS_solve_SCF_potmix_v2.jl")

function test_main()

    #Ham = create_Ham_H2()
    #Ham = create_Ham_Si_fcc()
    #Ham = create_Ham_Si_fcc(xcfunc="PBE")
    #Ham = create_Ham_GaAs_v1()
    #Ham = create_Ham_GaAs_v2()
    #Ham = create_Ham_CO()
    
    #Ham = create_Ham_N2()
    Ham = create_Ham_NH3()

    Random.seed!(1234)
    @time KS_solve_SCF_rhomix_v2!( Ham, mix_method="linear_adaptive", betamix=0.1 )

    #Random.seed!(1234)
    #@time KS_solve_SCF_rhomix_v2!( Ham, mix_method="simple", betamix=0.1 )

    #Random.seed!(1234)
    #@time KS_solve_SCF_rhomix_v2!(Ham, mix_method="broyden", betamix=0.5, mixdim=8)

    #Random.seed!(1234)
    #@time KS_solve_SCF_potmix!(Ham, mix_method="broyden", betamix=0.5, etot_conv_thr=1e-6)

    #Random.seed!(1234)
    #@time KS_solve_SCF!(Ham, mix_method="broyden", betamix=0.5, etot_conv_thr=1e-6)

    #Random.seed!(1234)
    #@time KS_solve_SCF!(Ham, mix_method="anderson", betamix=0.5, etot_conv_thr=1e-6)

    #Random.seed!(1234)
    #@time KS_solve_SCF!(Ham, mix_method="pulay", betamix=0.5, etot_conv_thr=1e-6)

    #Random.seed!(1234)
    #@time KS_solve_SCF!(Ham, mix_method="ppulay", betamix=0.5, etot_conv_thr=1e-6)

    #Random.seed!(1234)
    #@time KS_solve_Emin_PCG!(Ham)

end


test_main()
