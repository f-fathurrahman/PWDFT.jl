using Printf
using LinearAlgebra

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("create_Ham.jl")
include("symmetry_atoms.jl")

function main()
    Ham = create_Ham_Si_fcc()
    #Ham = create_Ham_GaAs_v1()

    init_irt( Ham.atoms, Ham.sym_info )
end

main()