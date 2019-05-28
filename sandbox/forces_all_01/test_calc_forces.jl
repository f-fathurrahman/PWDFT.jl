using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("../forces_ewald_01/calc_forces_NN.jl")
include("../forces_local_01/calc_forces_Ps_loc.jl")
include("../forces_nonlocal_01/calc_forces_Ps_nloc.jl")

function create_Ham_Si_fcc()

    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(5.431*ANG2BOHR))

    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]

    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[8,8,8] )
end



function create_Ham_GaAs_v1()

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
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[8,8,8] )
end


function create_Ham_GaAs_v2()
    atoms = Atoms(xyz_string_frac=
        """
        2

        Ga  0.0  0.0  0.0
        As  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(5.6537*ANG2BOHR))

    pspfiles = [joinpath(DIR_PSP, "Ga-q3.gth"),
                joinpath(DIR_PSP, "As-q5.gth")]

    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[8,8,8] )
end


function test_main()

    Ham = create_Ham_Si_fcc()
    #Ham = create_Ham_GaAs_v1()
    #Ham = create_Ham_GaAs_v2()

    println(Ham)

    Random.seed!(1234)
    
    KS_solve_Emin_PCG!(Ham, etot_conv_thr=1e-8, savewfc=true)
    #KS_solve_SCF!(Ham, mix_method="rpulay", etot_conv_thr=1e-8, savewfc=true)

    psiks = read_psiks(Ham)

    for psi in psiks
        ortho_check(psi)
    end

    atoms = Ham.atoms
    Natoms = atoms.Natoms
    atsymbs = atoms.atsymbs

    F_NN = calc_forces_NN( Ham.atoms )*2.0
    println("NN forces:")
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_NN[1,ia], F_NN[2,ia], F_NN[3,ia] )
    end

    F_Ps_loc = calc_forces_Ps_loc( Ham.atoms, Ham.pw, Ham.pspots, Ham.rhoe )*2.0
    println("Ps loc forces:")
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_Ps_loc[1,ia], F_Ps_loc[2,ia], F_Ps_loc[3,ia] )
    end

    F_Ps_nloc =
    calc_forces_Ps_nloc( Ham.atoms, Ham.pw, Ham.pspots, Ham.electrons, Ham.pspotNL, psiks )*2.0
    println("Ps nloc forces:")
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_Ps_nloc[1,ia], F_Ps_nloc[2,ia], F_Ps_nloc[3,ia] )
    end


    F_total = F_NN + F_Ps_loc + F_Ps_nloc
    println("Total forces:")
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_total[1,ia], F_total[2,ia], F_total[3,ia] )
    end

    @printf("Sum of forces in x-dir: %18.10f\n", sum(F_total[1,:]))
    @printf("Sum of forces in y-dir: %18.10f\n", sum(F_total[2,:]))
    @printf("Sum of forces in z-dir: %18.10f\n", sum(F_total[3,:]))
end


test_main()
