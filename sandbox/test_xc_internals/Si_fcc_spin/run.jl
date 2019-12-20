using LinearAlgebra
using Random
using Printf

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

import InteractiveUtils
InteractiveUtils.versioninfo()

function main()

    Random.seed!(1234)

    # Atoms
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3],
        use_xc_internal=true, Nspin=2 )
    println(Ham)

    #
    # Solve the KS problem
    #
    KS_solve_SCF!( Ham, mix_method="rpulay", betamix=0.5 )
    
end

main()
@time main()
