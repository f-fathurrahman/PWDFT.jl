using Printf
using PWDFT

using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include(joinpath(DIR_PWDFT, "sandbox", "ABINIT.jl"))
include(joinpath(DIR_PWDFT, "sandbox", "PWSCF.jl"))
include(joinpath(DIR_PWDFT, "sandbox", "KS_solve_SCF_potmix.jl"))

function init_Ham_GaAs_fcc( a::Float64, meshk::Array{Int64,1} )

    # Atoms
    atoms = Atoms( xyz_string_frac=
        """
        2

        Ga  0.0   0.0   0.0
        As  0.25  0.25  0.25
        """, in_bohr=true, LatVecs = gen_lattice_fcc(a) )

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Ga-q3.gth"),
                joinpath(DIR_PSP, "As-q5.gth")]
    ecutwfc = 15.0
    return Hamiltonian( atoms, pspfiles, ecutwfc, meshk=meshk )
end


function main()

    Random.seed!(1234)

    LATCONST = 10.6839444516

    scal_i = 0.8
    scal_f = 1.2
    Ndata = 15
    Δ = (scal_f - scal_i)/(Ndata-1)
    
    stdout_orig = stdout

    f = open("TEMP_EOS_data.dat", "w")

    for i = 1:Ndata

        scal = scal_i + (i-1)*Δ

        a = LATCONST*scal

        fname = "TEMP_LOG_"*string(i)
        fileout = open(fname, "w")
        
        redirect_stdout(fileout)
    
        Ham = init_Ham_GaAs_fcc( a, [3,3,3] )
        KS_solve_SCF!( Ham, mix_method="rpulay" )

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

        close(fileout)

        @printf(f, "%18.10f %18.10f %18.10f %18.10f\n", a, sum(Ham.energies),
                sum(abinit_energies), sum(pwscf_energies) )
    end

    close(f)

end

main()

