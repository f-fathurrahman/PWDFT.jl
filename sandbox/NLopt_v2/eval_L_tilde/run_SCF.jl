using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("../create_Ham.jl")

function main()

    #Ham = create_Ham_atom_Pt_smearing(a=10.0)
    #Ham = create_Ham_atom_Pt()
    Ham = create_Ham_Al_fcc_smearing(meshk=[1,1,1], Nspin=2)
    #Ham = create_Ham_Pt_fcc_smearing(meshk=[1,1,1])

    println(Ham)

    KS_solve_SCF!( Ham,
        print_final_ebands=true,
        mix_method="rpulay",
        betamix=0.2,
        use_smearing=true,
        savewfc=true )

    #KS_solve_Emin_PCG!( Ham, print_final_ebands=true, i_cg_beta=4 )

end

main()
