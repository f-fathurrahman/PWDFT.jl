using Printf
using LinearAlgebra
using Random

using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("KS_solve_SCF_04.jl")

include(joinpath(DIR_PWDFT, "sandbox", "ABINIT.jl"))
include(joinpath(DIR_PWDFT, "sandbox", "PWSCF.jl"))

function init_Ham_Si_fcc()
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0   0.0   0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true,
        LatVecs=gen_lattice_fcc(5.431*ANG2BOHR))

    ecutwfc = 15.0
    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
end


function init_Ham_GaAs()
    atoms = Atoms(xyz_string_frac=
        """
        2

        Ga  0.0  0.0  0.0
        As  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.6839444516))

    pspfiles = [joinpath(DIR_PSP, "Ga-q3.gth"),
                joinpath(DIR_PSP, "As-q5.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
end


function init_Ham_atom_H()
    atoms = Atoms( xyz_string="""
            1

            H   0.1   0.0   0.0
            """, LatVecs=gen_lattice_sc(16.0) )

    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc )
end

function init_Ham_H2()
    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "H2.xyz"),
                   LatVecs = gen_lattice_sc(16.0) )

    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc )
end

function init_Ham_H2_v2()
    atoms = Atoms( xyz_string=
        """
        2
            
        H  0.0  0.0  0.0
        H  1.5  0.0  0.0
        """, in_bohr=true, LatVecs = gen_lattice_sc(16.0) )

    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc )
end

function main()

    Random.seed!(1234)

    #Ham = init_Ham_Si_fcc()
    #Ham = init_Ham_GaAs()
    #Ham = init_Ham_atom_H()
    Ham = init_Ham_H2()
    #Ham = init_Ham_H2_v2()
    println(Ham)
    println(Ham.sym_info)

    @time KS_solve_SCF_04!( Ham, mix_method="anderson" )
    @time KS_solve_SCF_04!( Ham, mix_method="anderson" )
   
    run(`rm -fv TEMP_abinit/\*`)
    write_abinit(Ham, prefix_dir="./TEMP_abinit/")
    cd("./TEMP_abinit")
    run(pipeline(`abinit`, stdin="FILES", stdout="ABINIT_o_LOG"))
    cd("../")

    abinit_energies = read_abinit_etotal("TEMP_abinit/LOG1")
    println("\nABINIT result\n")
    println(abinit_energies)

    run(`rm -rfv TEMP_pwscf/\*`)
    write_pwscf( Ham, prefix_dir="TEMP_pwscf" )
    cd("./TEMP_pwscf")
    run(pipeline(`pw.x`, stdin="PWINPUT", stdout="LOG1"))
    cd("../")

    pwscf_energies = read_pwscf_etotal("TEMP_pwscf/LOG1")
    println("\nPWSCF result\n")
    println(pwscf_energies)

end

main()
