using Printf
using LinearAlgebra
using Random
using PWDFT
import Serialization

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT,"pseudopotentials","pade_gth")

include("../common/dump_bandstructure.jl")
include("../common/gen_kpath.jl")

function main()

    Random.seed!(1234)

    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true,
        LatVecs=gen_lattice_fcc(10.2631))
    
    pspfiles = [joinpath(DIR_PSP,"Si-q4.gth")]
    ecutwfc = 15.0

    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], extra_states=4 )
    KS_solve_SCF!( Ham, mix_method="rpulay" )

    Serialization.serialize("Ham.data", Ham)

end

main()
