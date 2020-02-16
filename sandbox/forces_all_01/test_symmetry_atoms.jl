using Printf
using LinearAlgebra

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("create_Ham.jl")
include("symmetry_atoms.jl")

function main()
    #Ham = create_Ham_CO()
    Ham = create_Ham_Si_fcc()
    #Ham = create_Ham_GaAs_v1()

    irt = init_irt( Ham.atoms, Ham.sym_info )

    v = zeros(3,2)
    #v[:,1] = [1.0, 1.0, 1.0]
    #v[:,2] = [1.0, 0.0, 0.0]
    v[:,1] = [-0.1, 1.0, 0.0]
    v[:,2] = [ 0.0, 0.10, 1.0]


    display(v); println()

    symmetrize_vector!( Ham.pw, Ham.sym_info, irt, v)

    display(v); println()

end

main()