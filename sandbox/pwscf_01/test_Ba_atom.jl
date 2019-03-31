using Printf
using PWDFT

const PWDFT_DIR = joinpath(dirname(pathof(PWDFT)), "..")
const PSPOT_PATH = joinpath(PWDFT_DIR, "pseudopotentials", "pade_gth")

include(joinpath(PWDFT_DIR, "sandbox", "PWSCF.jl"))

function test_main()

    atoms = Atoms( xyz_string="""
        1

        Ba  0.0  0.0  0.0
        """, LatVecs=gen_lattice_sc(16.0))

    pspfiles = [joinpath(PSPOT_PATH, "Ba-q10.gth")]
    
    ecutwfc = 15.0
    Ham = Hamiltonian(atoms, pspfiles, ecutwfc, extra_states=4)

    #KS_solve_SCF!( Ham, use_smearing=true, kT=0.01, mix_method="rpulay" )
    
    write_pwscf( Ham, use_smearing=true, kT=0.01, prefix_dir="TEMP_pwscf" )
    
    run(`rm -rfv TEMP_pwscf/\*`)  # need recursive option for rm
    write_pwscf( Ham, use_smearing=true, kT=0.01, prefix_dir="TEMP_pwscf" )
    cd("./TEMP_pwscf")
    run(pipeline(`pw.x`, stdin="PWINPUT", stdout="LOG1"))
    cd("../")

    energies = read_pwscf_etotal("TEMP_pwscf/LOG1")
    println(energies)

end

test_main()
