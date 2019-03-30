using Printf
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include(joinpath(DIR_PWDFT, "sandbox", "ABINIT.jl"))
include(joinpath(DIR_PWDFT, "sandbox", "KS_solve_SCF_potmix.jl"))


function init_Ham_Al_fcc( a::Float64, meshk::Array{Int64,1} )
    # Atoms
    atoms = Atoms( xyz_string_frac=
        """
        1

        Al  0.0  0.0  0.0
        """, in_bohr=true,
        LatVecs = gen_lattice_fcc(a) )
    pspfiles = [joinpath(DIR_PSP, "Al-q3.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, xcfunc="LDA",
                       Nspin=1, meshk=[8,8,8], extra_states=4 )
end




function main()

    LATCONST = 7.6525970200

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
    
        Ham = init_Ham_Al_fcc( a, [8,8,8] )
        KS_solve_SCF!( Ham, use_smearing=true, kT=0.001, mix_method="rpulay" )

        run(`rm -fv TEMP_abinit/\*`)
        write_abinit(Ham, prefix_dir="./TEMP_abinit/", use_smearing=true, kT=0.001)
        cd("./TEMP_abinit")
        run(pipeline(`abinit`, stdin="FILES", stdout="ABINIT_o_LOG"))
        cd("../")

        abinit_energies = read_abinit_etotal("TEMP_abinit/LOG1")
        println("\nABINIT result\n")
        println(abinit_energies)

        close(fileout)

        @printf(f, "%18.10f %18.10f %18.10f\n", a, sum(Ham.energies), sum(abinit_energies) )
    end

    close(f)

end

main()

