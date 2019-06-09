using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("precKerker.jl")
include("mix_broyden.jl")
include("create_Ham.jl")
include("KS_solve_SCF_rhomix_v2.jl")
include("KS_solve_SCF_potmix_v2.jl")

function test_main()

    #Ham = create_Ham_O2_smearing()
    Ham = create_Ham_Pt_fcc_smearing()

    Random.seed!(1234)
    @time KS_solve_SCF_rhomix_v2!(Ham, mix_method="rpulay", betamix=0.5, mixdim=5, use_smearing=true)
    Random.seed!(1234)
    @time KS_solve_SCF_rhomix_v2!(Ham, mix_method="rpulay", betamix=0.5, mixdim=5, use_smearing=true)

    Random.seed!(1234)
    @time KS_solve_SCF_potmix_v2!(Ham, mix_method="simple", betamix=0.5, mixdim=8, use_smearing=true)
    Random.seed!(1234)
    @time KS_solve_SCF_potmix_v2!(Ham, mix_method="simple", betamix=0.5, mixdim=8, use_smearing=true)

end


test_main()
