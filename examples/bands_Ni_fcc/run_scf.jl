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
    
    atoms = Atoms( xyz_string_frac=
        """
        1

        Ni  0.0  0.0  0.0
        """, in_bohr=true,
        LatVecs = gen_lattice_fcc(6.65914911201) )

    pspots = [PsPot_GTH(joinpath(DIR_PSP,"Ni-q10.gth"))]
    ecutwfc = 50.0

    options = HamiltonianOptions()
    options.meshk = [8,8,8]
    options.extra_states = 4
    options.Nspin_channel = 2
    options.Nspin_comp = 2
    Ham = Hamiltonian( atoms, pspots, ecutwfc, options )
    KS_solve_SCF!( Ham, mix_method="anderson",
        betamix=0.1, starting_magn=[0.1], use_smearing=true
    )

    Serialization.serialize("Ham_ecut_"*string(ecutwfc)*".data", Ham)

end

main()
