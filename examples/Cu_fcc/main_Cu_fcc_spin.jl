function main()
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        1

        Cu  0.0  0.0  0.0
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(3.61496*ANG2BOHR)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)

    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/Cu-q11.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, Nspin=2, xcfunc="VWN",
                         meshk=[3,3,3], verbose=true, extra_states=4 )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

    #
    # Solve the KS problem
    #
    KS_solve_SCF_smearing!( Ham, β=0.5, mix_method="anderson" )

    println("\nTotal energy components")
    println(Ham.energies)

end
