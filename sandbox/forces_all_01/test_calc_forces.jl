using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("calc_forces_NN.jl")
include("calc_forces_Ps_loc.jl")
include("calc_forces_Ps_nloc.jl")

#include("symmetry_atoms.jl")
#include("../symmetry_qe/SymmetryBase.jl")

include("create_Ham.jl")


function test_main()

    #Ham = create_Ham_H2()
    Ham = create_Ham_Si_fcc()
    #Ham = create_Ham_Si_fcc(xcfunc="PBE")
    #Ham = create_Ham_GaAs_v1()
    #Ham = create_Ham_GaAs_v2()
    #Ham = create_Ham_CO()
    println(Ham)

    #sym_base = SymmetryBase(Ham.atoms)
    #println(Ham.sym_info)
    #irt = init_irt(Ham.atoms, Ham.sym_info)

    Random.seed!(1234)
    
    #KS_solve_Emin_PCG!(Ham, savewfc=true)
    #psiks = read_psiks(Ham)

    psiks = rand_BlochWavefunc(Ham)
    KS_solve_SCF!( Ham, psiks, mix_method="anderson", etot_conv_thr=1e-6 )

    atoms = Ham.atoms
    Natoms = atoms.Natoms
    atsymbs = atoms.atsymbs

    # All forces are multiplied by 2 to match QE unit (Ry/bohr)

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
    
    PWDFT.symmetrize_vector!( Ham.pw.LatVecs, Ham.sym_info, F_Ps_nloc )
    println("Ps nloc forces symmetrized:")
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_Ps_nloc[1,ia], F_Ps_nloc[2,ia], F_Ps_nloc[3,ia] )
    end



    F_total = F_NN + F_Ps_loc + F_Ps_nloc
    println("Total forces: (in Ry unit)")
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_total[1,ia], F_total[2,ia], F_total[3,ia] )
    end

    @printf("Sum of forces in x-dir: %18.10f\n", sum(F_total[1,:]))
    @printf("Sum of forces in y-dir: %18.10f\n", sum(F_total[2,:]))
    @printf("Sum of forces in z-dir: %18.10f\n", sum(F_total[3,:]))

    F_total = 0.5*F_total
    println("Total forces (Hartree unit):")
    for ia = 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
                F_total[1,ia], F_total[2,ia], F_total[3,ia] )
    end


    #symmetrize_vector!( Ham.pw, Ham.sym_info, irt, F_total )
    #println("Total forces symmetrized:")
    #for ia = 1:Natoms
    #    @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia],
    #            F_total[1,ia], F_total[2,ia], F_total[3,ia] )
    #end

end


test_main()
