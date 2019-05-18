using Printf
using PWDFT
using LinearAlgebra

using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("rand_wfc.jl")

function init_Ham_Si_fcc()
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))

    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
end

function main()
    Random.seed!(1234)
    Ham = init_Ham_Si_fcc()

    psiks = rand_wfc(Ham)
end

main()

