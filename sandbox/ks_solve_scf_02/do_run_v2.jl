using Printf
using PWDFT

using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include(joinpath(DIR_PWDFT, "sandbox", "ABINIT.jl"))
include(joinpath(DIR_PWDFT, "sandbox", "PWSCF.jl"))
include(joinpath(DIR_PWDFT, "sandbox", "KS_solve_SCF_potmix.jl"))
include("KS_solve_SCF_02.jl")

function init_Ham_Si_fcc( a::Float64, meshk::Array{Int64,1} )
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(a))

    pspfiles = ["Si-q4-local-only.gth"]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=meshk, Ns_=(32,32,32) )
end

function main()

    Random.seed!(1234)

    LATCONST = 10.2631

    Ham = init_Ham_Si_fcc( LATCONST, [8,8,8] )
    println(Ham)
    println("sum V_ps_loc = ", sum(Ham.potentials.Ps_loc)*2)  # convert to Ry
    println("V ps loc = ", Ham.potentials.Ps_loc[1:5]*2)
    #exit()

    #KS_solve_Emin_PCG!( Ham )
    KS_solve_SCF_02!( Ham, betamix=0.1 )

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

