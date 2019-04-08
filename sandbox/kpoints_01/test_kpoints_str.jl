using Random
using LinearAlgebra

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

function main()

    LATCONST = 10.2631

    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0   0.0   0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(LATCONST))

    ecutwfc = 15.0
    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    kpts_str = """
    4
    0.000000000000   0.00000000000   0.00000000000   0.074074074074
    0.000000000000   0.00000000000   0.33333333333   0.592592592592
    0.000000000000   0.33333333333   0.33333333333   0.444444444444
    0.000000000000   0.33333333333  -0.33333333333   0.888888888889
    """

    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, Ns_=(32,32,32), kpts_str=kpts_str )

    println(Ham)

    #KS_solve_SCF!( Ham, mix_method="rpulay" )
    KS_solve_Emin_PCG!( Ham )

end

main()

