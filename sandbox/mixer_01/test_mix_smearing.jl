using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("precKerker.jl")
include("mix_adaptive.jl")
include("create_Ham.jl")
include("KS_solve_SCF_rhomix_v2.jl")
include("KS_solve_SCF_potmix_v2.jl")

function test_main()

    #Ham = create_Ham_O2_smearing()
    #Ham = create_Ham_Pt_fcc_smearing()
    Ham = create_Ham_Fe_bcc()
    #Ham = create_Ham_Si_fcc( Nspin=2 )

    @time KS_solve_SCF_rhomix_v2!( Ham, mix_method="linear_adaptive",
        betamix=0.1, use_smearing=true, starting_magnetization=[0.5] )

    @time KS_solve_SCF_rhomix_v2!( Ham, mix_method="simple",
        betamix=0.1, use_smearing=true, starting_magnetization=[0.5] )

    #Random.seed!(1234)
    #@time KS_solve_SCF_rhomix_v2!(Ham, mix_method="rpulay", betamix=0.5, mixdim=5, use_smearing=true)
    #Random.seed!(1234)
    #@time KS_solve_SCF_rhomix_v2!(Ham, mix_method="rpulay", betamix=0.5, mixdim=5, use_smearing=true)

    #Random.seed!(1234)
    #@time KS_solve_SCF_potmix_v2!(Ham, mix_method="simple", betamix=0.5, mixdim=8, use_smearing=true)
    #Random.seed!(1234)
    #@time KS_solve_SCF_potmix_v2!(Ham, mix_method="simple", betamix=0.5, mixdim=8, use_smearing=true)

    #Random.seed!(1234)
    #@time KS_solve_SCF_potmix!(Ham, mix_method="simple", betamix=0.5, use_smearing=true)
    
    #Random.seed!(1234)
    #@time KS_solve_SCF_potmix!(Ham, mix_method="simple", betamix=0.5, use_smearing=true, starting_magnetization=[0.5])

    #Random.seed!(1234)
    #@time KS_solve_SCF_potmix!(Ham, mix_method="broyden", betamix=0.5, use_smearing=true, starting_magnetization=[0.5])

    #Random.seed!(1234)
    #@time KS_solve_SCF!(Ham, mix_method="broyden", betamix=0.1, use_smearing=true)

    #Random.seed!(1234)
    #@time KS_solve_SCF!(Ham, mix_method="anderson", betamix=0.5, use_smearing=true, starting_magnetization=[0.5])

    #Random.seed!(1234)
    #@time KS_solve_SCF!(Ham, mix_method="pulay", betamix=0.5, use_smearing=true, starting_magnetization=[0.5])
    
    #Random.seed!(1234)
    #@time KS_solve_SCF!(Ham, mix_method="rpulay", betamix=0.5, use_smearing=true, starting_magnetization=[0.5])

    #Random.seed!(1234)
    #@time KS_solve_SCF!(Ham, mix_method="ppulay", betamix=0.5, use_smearing=true, starting_magnetization=[0.5])

    #Random.seed!(1234)
    #@time KS_solve_SCF!(Ham, mix_method="broyden", betamix=0.5, use_smearing=true, starting_magnetization=[0.5])

end


test_main()
