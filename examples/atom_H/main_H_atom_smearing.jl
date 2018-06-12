using PWDFT

function test_main()
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        1

        H  0.0  0.0  0.0
        """)
    atoms.LatVecs = gen_lattice_cubic(16.0)
    println(atoms)

    # Initialize Hamiltonian
    #pspfiles = ["../../pseudopotentials/pade_gth/H-q1.gth"]
    pspfiles = ["../../pseudopotentials/pbe_gth/H-q1.gth"]
    ecutwfc_Ry = 30.0
    Ham = PWHamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, xcfunc="PBE", Nspin=2, extra_states=4 )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

    KS_solve_SCF_smearing!( Ham, mix_method="simple",  update_psi="PCG", β=0.1, ETOT_CONV_THR=1e-6 )

    Nstates = Ham.electrons.Nstates
    ebands = Ham.electrons.ebands
    
    println("\nBand energies:")
    for ist = 1:Nstates
        @printf("%8d  %18.10f = %18.10f eV\n", ist, ebands[ist], ebands[ist]*Ry2eV*2)
    end
    
    println("\nTotal energy components")
    println(Ham.energies)

end

@time test_main()
