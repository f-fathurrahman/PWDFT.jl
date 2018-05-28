using PWDFT

function test_main()

    # Atoms
    atoms = init_atoms_xyz_string(
    """
    1
    
    Ni   0.0   0.0   0.0
    """)
 
    atoms.LatVecs = gen_lattice_cubic(16.0)
    println(atoms)

    #
    # Initialize Hamiltonian
    #
    pspfiles = ["../../pseudopotentials/pade_gth/Ni-q18.gth"]
    ecutwfc_Ry = 30.0
    Ham = PWHamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, verbose=true,
                         Nspin=1 )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

    #KS_solve_SCF!( Ham, mix_method="anderson", β=0.2 )
    KS_solve_Emin_PCG!( Ham )

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

