using LinearAlgebra
using Printf
using Random

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include(joinpath(DIR_PWDFT, "sandbox", "KS_solve_SCF_NLsolve.jl"))

include("calc_energies_grad.jl")
include("create_Ham.jl")

function main()
    Random.seed!(1234)

    #Ham = create_Ham_H2()
    #Ham = create_Ham_H_atom()
    #Ham = create_Ham_Si_fcc()
    Ham = create_Ham_GaAs()
    #Ham = create_Ham_NH3()
    #Ham = create_Ham_ZnO()

    psiks = rand_BlochWavefunc( Ham )

    Rhoe = calc_rhoe(Ham, psiks)
    update!(Ham, Rhoe)

    g = zeros_BlochWavefunc( Ham )
    Kg = zeros_BlochWavefunc( Ham )

    Ham.energies.NN = calc_E_NN( Ham.atoms )
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )

    E0 = calc_energies_grad!( Ham, psiks, g, Kg )
    println("E0 = ", E0)

    psic = zeros_BlochWavefunc( Ham )
    for Δ in 10.0 .^ range(1,stop=-9,step=-1)
        psic = deepcopy( psiks )
        do_step!( psiks, )
    end

end

main()