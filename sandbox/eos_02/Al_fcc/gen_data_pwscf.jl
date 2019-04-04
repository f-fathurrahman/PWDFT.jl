using Printf
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include(joinpath(DIR_PWDFT, "sandbox", "ABINIT.jl"))
include(joinpath(DIR_PWDFT, "sandbox", "PWSCF.jl"))
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
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, xcfunc="VWN",
                       Nspin=1, meshk=meshk, extra_states=4 )
end



function main()

    LATCONST = 7.6525970200

    scal_i = 0.8
    scal_f = 1.2
    Ndata = 15
    Δ = (scal_f - scal_i)/(Ndata-1)

    f = open("TEMP_EOS_data_pwscf.dat", "w")

    for i = 1:Ndata

        scal = scal_i + (i-1)*Δ

        a = LATCONST*scal

        Ham = init_Ham_Al_fcc( a, [8,8,8] )

        run(`rm -rfv TEMP_pwscf/\*`)
        write_pwscf( Ham, prefix_dir="TEMP_pwscf", use_smearing=true, kT=0.001 )
        cd("./TEMP_pwscf")
        run(pipeline(`pw.x`, stdin="PWINPUT", stdout="LOG1"))
        cd("../")

        pwscf_energies = read_pwscf_etotal("TEMP_pwscf/LOG1")
        println("\nPWSCF result\n")
        println(pwscf_energies)

        @printf(f, "%18.10f %18.10f\n", a, sum(pwscf_energies) )

    end

    close(f)

end

main()

