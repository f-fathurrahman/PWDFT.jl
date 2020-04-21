using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("smearing.jl")
include("create_Ham.jl")
include("ElecVars.jl")
include("test_ElecVars.jl")
include("calc_energies_grad.jl")

function main()

    Random.seed!(1234)

    kT = 0.01
    #Ham = create_Ham_atom_Al_smearing()
    Ham = create_Ham_Al_fcc_smearing()
    #Ham = create_Ham_atom_Pt_smearing()
    #Ham = create_Ham_Pt_fcc_smearing()
    println(Ham)

    #test_ElecVars(Ham)

    evars = ElecVars(Ham)
    g = ElecGradient(Ham)
    Kg = ElecGradient(Ham)
    Ham.energies.NN = calc_E_NN(Ham.atoms)
    Etot = calc_energies_grad!( Ham, evars, g, Kg, kT )

    println(Ham.energies)
    println("Etot = ", Etot)

    println("Pass here")
end

main()