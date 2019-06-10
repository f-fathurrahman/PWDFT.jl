using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("calc_forces_Ps_nloc.jl")

function test_Si_fcc()

    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.1  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(5.431*ANG2BOHR))

    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]

    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )

    println(Ham)

    Random.seed!(1234)
    KS_solve_Emin_PCG!( Ham, savewfc=true )
    #KS_solve_SCF!(Ham, mix_method="rpulay", etot_conv_thr=1e-6, savewfc=true)

    psiks = read_psiks(Ham)

    Natoms = atoms.Natoms
    atsymbs = atoms.atsymbs

#    println("")
#    F_Ps_loc = calc_forces_Ps_nloc( Ham.atoms, Ham.pw, Ham.pspots, Ham.rhoe )*2.0
#    for ia = 1:Natoms
#        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
#                F_Ps_loc[1,ia], F_Ps_loc[2,ia], F_Ps_loc[3,ia] )
#    end

    println("")
    F_Ps_nloc =
    calc_forces_Ps_nloc_finite_diff( Ham.atoms, Ham.pw, Ham.pspots, Ham.electrons, psiks )*2.0
    println("Using finite difference")
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_Ps_nloc[1,ia], F_Ps_nloc[2,ia], F_Ps_nloc[3,ia] )
    end

    #println("\nTest calc_dbetaNL:")
    #dbetaNL = calc_dbetaNL(Ham.atoms, Ham.pw, Ham.pspots, Ham.pspotNL)
    #println("sum dbetaNL = ", sum(dbetaNL))

    #ik = 1
    #dbetaNL_psi = calc_dbetaNL_psi(ik, dbetaNL, psiks[1])
    #println("sum dbetaNL_psi = ", sum(dbetaNL_psi))

    F_Ps_nloc =
    calc_forces_Ps_nloc( Ham.atoms, Ham.pw, Ham.pspots, Ham.electrons, Ham.pspotNL, psiks )*2.0
    println("Analytic:")
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_Ps_nloc[1,ia], F_Ps_nloc[2,ia], F_Ps_nloc[3,ia] )
    end


end

#=
function test_GaAs_v1()

    LatVecs = zeros(3,3)
    LatVecs[:,1] = [0.5, 0.5, 0.0]
    LatVecs[:,2] = [0.5, 0.0, 0.5]
    LatVecs[:,3] = [0.0, 0.5, 0.5]
    LatVecs = LatVecs*5.6537*ANG2BOHR

    atoms = Atoms(xyz_string_frac=
        """
        2

        Ga  0.0  0.0  0.0
        As  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=LatVecs)

    pspfiles = [joinpath(DIR_PSP, "Ga-q3.gth"),
                joinpath(DIR_PSP, "As-q5.gth")]

    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )

    Random.seed!(1234)
    #KS_solve_Emin_PCG!(Ham, etot_conv_thr=1e-8)
    KS_solve_SCF!(Ham, mix_method="rpulay", etot_conv_thr=1e-6)

    Natoms = atoms.Natoms
    atsymbs = atoms.atsymbs

#    println("")
#    F_Ps_loc = calc_forces_Ps_loc( Ham.atoms, Ham.pw, Ham.pspots, Ham.rhoe )*2.0
#    for ia = 1:Natoms
#        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
#                F_Ps_loc[1,ia], F_Ps_loc[2,ia], F_Ps_loc[3,ia] )
#    end

#    println("")
#    F_Ps_loc = calc_forces_Ps_loc_finite_diff( Ham.atoms, Ham.pw, Ham.pspots, Ham.rhoe )*2.0
#    println("Using finite difference")
#    for ia = 1:Natoms
#        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
#                F_Ps_loc[1,ia], F_Ps_loc[2,ia], F_Ps_loc[3,ia] )
#    end

end
=#


#=
function test_GaAs_v2()

    atoms = Atoms(xyz_string_frac=
        """
        2

        Ga  0.0  0.0  0.0
        As  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(5.6537*ANG2BOHR))

    pspfiles = [joinpath(DIR_PSP, "Ga-q3.gth"),
                joinpath(DIR_PSP, "As-q5.gth")]

    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )

    Random.seed!(1234)
    #KS_solve_Emin_PCG!( Ham, etot_conv_thr=1e-8 )
    KS_solve_SCF!( Ham, mix_method="rpulay", etot_conv_thr=1e-6, savewfc=true )

    psiks = read_psiks( Ham )

    Natoms = atoms.Natoms
    atsymbs = atoms.atsymbs

    println("")
    F_Ps_loc = calc_forces_Ps_nloc( Ham.atoms, Ham.pw, Ham.pspots, Ham.rhoe )*2.0
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_Ps_loc[1,ia], F_Ps_loc[2,ia], F_Ps_loc[3,ia] )
    end

    println("")
    F_Ps_loc = calc_forces_Ps_loc_finite_diff( Ham.atoms, Ham.pw, Ham.pspots, Ham.rhoe )*2.0
    println("Using finite difference")
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_Ps_loc[1,ia], F_Ps_loc[2,ia], F_Ps_loc[3,ia] )
    end



end
=#

#test_H2()
test_Si_fcc()
#test_GaAs_v1()
#test_GaAs_v2()
