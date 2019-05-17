using Printf
using Random
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("guess_rhoe_atomic.jl")

function init_Ham_GaAs()
    atoms = Atoms(xyz_string_frac=
        """
        2

        Ga  0.0  0.0  0.0
        As  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.6839444516))

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Ga-q3.gth"),
                joinpath(DIR_PSP, "As-q5.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], Nspin=2 )
end

function main()
    Ham = init_Ham_GaAs()
    Rhoe = guess_rhoe_atomic( Ham, starting_magnetization=[0.0, 0.1])

    dVol = Ham.pw.CellVolume/prod(Ham.pw.Ns)
    magn = Rhoe[:,1] - Rhoe[:,2]
    println("integ magn = ", sum(magn)*dVol)
    println("integ sum Rhoe = ", sum(Rhoe)*dVol)
end

main()
