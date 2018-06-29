using Printf
using PWDFT

function test_main()
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        1

        Cu  0.0  0.0  0.0
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(3.61496*ANG2BOHR)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)
    write_xsf( "TEMP_Cu.xsf", atoms )

    # Initialize Hamiltonian
    pspfiles = ["../../pseudopotentials/pade_gth/Cu-q11.gth"]
    ecutwfc_Ry = 30.0
    Ham = PWHamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5,
                         meshk=[8,8,8], verbose=true, extra_states=4 )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

    #
    # Solve the KS problem
    #
    KS_solve_SCF_smearing!( Ham, β=0.2, mix_method="anderson", NiterMax=50 )

    println("\nTotal energy components")
    println(Ham.energies)

end

@time test_main()


