using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("calc_forces_Ps_loc.jl")

function test_H2()

    atoms = Atoms(xyz_string=
        """
        2

        H      3.83653478       4.23341768       4.23341768
        H      4.63030059       4.23341768       4.23341768
        """, LatVecs=gen_lattice_sc(16.0))

    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]

    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc )

    Random.seed!(1234)
    KS_solve_Emin_PCG!(Ham, ETOT_CONV_THR=1e-8)

    Natoms = atoms.Natoms
    atsymbs = atoms.atsymbs

    println("")
    F_Ps_loc = calc_forces_Ps_loc( Ham )*2.0
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_Ps_loc[1,ia], F_Ps_loc[2,ia], F_Ps_loc[3,ia] )
    end

#    F_Ps_loc = calc_forces_Ps_loc_finite_diff( Ham )*2.0
#    println("")
#    println("Using finite difference")
#    for ia = 1:Natoms
#        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
#                F_Ps_loc[1,ia], F_Ps_loc[2,ia], F_Ps_loc[3,ia] )
#    end

end

test_H2()
