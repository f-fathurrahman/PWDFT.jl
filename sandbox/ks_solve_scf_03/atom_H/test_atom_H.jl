using Printf
using PWDFT
using LinearAlgebra

using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include(joinpath(DIR_PWDFT, "sandbox", "ABINIT.jl"))
include(joinpath(DIR_PWDFT, "sandbox", "PWSCF.jl"))

include("../KS_solve_SCF_03.jl")

function init_Ham_atom_H()
    # Atoms
    atoms = Atoms( xyz_string="""
            1

            H   0.0   0.0   0.0
            """, LatVecs = gen_lattice_sc(16.0) )

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc )
end

function main()

    Random.seed!(1234)

    Ham = init_Ham_atom_H()


    #println("")
    #println("======================")
    #println("Using KS_solve_SCF:")
    #println("======================")
    #println("")
    #KS_solve_SCF!( Ham, mix_method="rpulay" )
    #e1 = copy( Ham.energies )

    println("")
    println("======================")
    println("Using KS_solve_SCF_03:")
    println("======================")
    println("")
    KS_solve_SCF_03!( Ham, mix_method="rpulay" )
    e2 = copy( Ham.energies )

    exit()


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


    println()
    println("Difference between KS_solve_SCF and KS_solve_SCF_03:")
    println()
    println( e1 - e2 )


    println()
    println("Difference between KS_solve_SCF and ABINIT:")
    println()
    println( e1 - abinit_energies )


    println()
    println("Difference between KS_solve_SCF_03 and ABINIT:")
    println()
    println( e2 - abinit_energies )

end

main()

