using Printf
using LinearAlgebra
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("create_Ham.jl")

include("KS_solve_DCM_01.jl")
include("KS_solve_TRDCM_01.jl")

function main()

    Ham = create_Ham_Si_fcc()
    #Ham = create_Ham_CO()
    #Ham = create_Ham_H2()

    #KS_solve_DCM_01!( Ham, NiterMax=20, MaxInnerSCF=3 )
    KS_solve_TRDCM_01!( Ham, NiterMax=20 )
    #KS_solve_Emin_PCG!( Ham )

    #psiks_old = rand_BlochWavefunc( Ham )
    #Rhoe = calc_rhoe( Ham, psiks_old )
#
#    #update!( Ham, Rhoe )
#
#    #psiks = rand_BlochWavefunc( Ham )
#    #@time test_shifted_Hamiltonian( 0.5, Ham, psiks_old, psiks )
    #@time test_shifted_Hamiltonian( 0.5, Ham, psiks_old, psiks )

    #my_trscf!( Ham )
end

main()