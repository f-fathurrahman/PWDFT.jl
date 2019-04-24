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

function init_Ham_Pt_fcc( a::Float64, meshk::Array{Int64,1} )
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, in_bohr=true, LatVecs=gen_lattice_fcc(a))

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Pt-q10.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=meshk, extra_states=4)
end

function main()

    Random.seed!(1234)

    LATCONST = 3.9231*ANG2BOHR

    Ham = init_Ham_Pt_fcc( LATCONST, [3,3,3] )
    println(Ham)

    println("")
    println("======================")
    println("Using KS_solve_SCF:")
    println("======================")
    println("")
    KS_solve_SCF!( Ham, mix_method="rpulay", use_smearing=true )
    e1 = copy( Ham.energies )

    println("")
    println("======================")
    println("Using KS_solve_SCF_03:")
    println("======================")
    println("")
    KS_solve_SCF_03!( Ham, mix_method="rpulay", use_smearing=true )
    e2 = copy( Ham.energies )


    run(`rm -fv TEMP_abinit/\*`)
    write_abinit(Ham, prefix_dir="./TEMP_abinit/", use_smearing=true )
    cd("./TEMP_abinit")
    run(pipeline(`abinit`, stdin="FILES", stdout="ABINIT_o_LOG"))
    cd("../")

    abinit_energies = read_abinit_etotal("TEMP_abinit/LOG1")
    println("\nABINIT result\n")
    println(abinit_energies)

    run(`rm -rfv TEMP_pwscf/\*`)
    write_pwscf( Ham, prefix_dir="TEMP_pwscf", use_smearing=true )
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

