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

    # Atoms
    atoms = Atoms( xyz_string_frac=
        """
        1

        Al  0.0  0.0  0.0
        """, in_bohr=true,
        LatVecs=gen_lattice_fcc(4.0495*ANG2BOHR) )

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP,"Al-q3.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc,
                       meshk=[8,8,8], extra_states=4 )
    println(Ham)

    # Solve the KS problem
    KS_solve_SCF!( Ham, betamix=0.2, mix_method="rpulay", use_smearing=true, kT=0.001 )

    Serialization.serialize("Ham.data", Ham)

end

main()
