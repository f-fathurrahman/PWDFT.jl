using Printf
using PWDFT

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)), "..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")
const DIR_STRUCTURES = joinpath(DIR_PWDFT, "structures")

include("../ABINIT.jl")

function test_CuSO4()
    # initialize atoms and Hamiltonian
    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "CuSO4.xyz"),
                   LatVecs=gen_lattice_sc(20.0) )
    pspfiles = [joinpath(DIR_PSP, "Cu-q11.gth"),
                joinpath(DIR_PSP, "O-q6.gth"),
                joinpath(DIR_PSP, "S-q6.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc )

    run(`rm -rfv TEMP_abinit/\*`)
    write_abinit( Ham, use_smearing=true, kT=0.01, prefix_dir="./TEMP_abinit/" )
    cd("./TEMP_abinit")
    run(pipeline(`abinit`, stdin="FILES", stdout="LOG1"))
    cd("../")

    energies = read_abinit_etotal("TEMP_abinit/LOG1")
    println(energies)
end


function test_H2()
    # initialize atoms and Hamiltonian
    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "H2.xyz"),
                   LatVecs=gen_lattice_sc(16.0) )
    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc )

    run(`rm -rfv TEMP_abinit/\*`)
    write_abinit( Ham, prefix_dir="./TEMP_abinit/" )
    cd("./TEMP_abinit")
    run(pipeline(`abinit`, stdin="FILES", stdout="LOG1"))
    cd("../")

    energies = read_abinit_etotal("TEMP_abinit/LOG1")
    println(energies)
end


function test_Al_fcc()

    atoms = Atoms( xyz_string_frac=
    """
    1

    Al  0.0  0.0  0.0
    """, in_bohr=true,
    LatVecs = gen_lattice_fcc(7.6525970200) )
    pspfiles = [joinpath(DIR_PSP, "Al-q3.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, xcfunc="LDA",
                       Nspin=1, meshk=[8,8,8], extra_states=4 )

    run(`rm -rfv TEMP_abinit/\*`)
    write_abinit( Ham, use_smearing=true, kT=0.01, prefix_dir="./TEMP_abinit/" )
    cd("./TEMP_abinit")
    run(pipeline(`abinit`, stdin="FILES", stdout="LOG1"))
    cd("../")

    energies = read_abinit_etotal("TEMP_abinit/LOG1")
    println(energies)

end


function test_Pt_fcc()
    # Atoms
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_fcc(3.9231*ANG2BOHR))

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Pt-q18.gth")]
    ecutwfc = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, xcfunc="PBE",
                        meshk=[3,3,3], extra_states=4 )
    
    run(`rm -rfv TEMP_abinit/\*`)
    write_abinit( Ham, use_smearing=true, kT=0.01, prefix_dir="TEMP_abinit/" )
    cd("./TEMP_abinit")
    run(pipeline(`abinit`, stdin="FILES", stdout="LOG1"))
    cd("../")

    energies = read_abinit_etotal("TEMP_abinit/LOG1")
    println(energies)

end

test_CuSO4()
test_H2()
test_Al_fcc()
test_Pt_fcc()
