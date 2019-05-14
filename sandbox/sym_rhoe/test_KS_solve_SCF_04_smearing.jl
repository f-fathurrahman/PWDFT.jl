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

function init_Ham_Pt_fcc()
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, in_bohr=true, LatVecs=gen_lattice_fcc(3.9231*ANG2BOHR))

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Pt-q10.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], extra_states=4)
end

function main()

    Random.seed!(1234)

    Ham = init_Ham_Pt_fcc()
    println(Ham)
    println(Ham.sym_info)

    @time KS_solve_SCF_04!( Ham, mix_method="rpulay", use_smearing=true )
    @time KS_solve_SCF_04!( Ham, mix_method="rpulay", use_smearing=true )
   
    run(`rm -fv TEMP_abinit/\*`)
    write_abinit( Ham, prefix_dir="./TEMP_abinit/", use_smearing=true )
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

end

main()
