using PWDFT

function test_main()
    const LatVecs = 16.0*diagm( ones(3) )
    ecutwfc_Ry = 40.0*0.5
    pw = PWGrid( ecutwfc_Ry, LatVecs )
    println(pw)

    #
    # Atoms
    #
    atoms = init_atoms_xyz("LiH.xyz")
    println(atoms)

    #
    # Structure factor
    #
    strf = calc_strfact( atoms, pw )

    #
    # Initialize Hamiltonian
    #
    Ham = PWHamiltonian(pw)
    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth",
                "../pseudopotentials/pade_gth/Li-q3.gth"]
    update!(Ham, atoms, strf, pspfiles)
    println("sum V Ps loc = ", sum(Ham.potentials.Ps_loc))

    #
    # calculate E_NN
    #
    Zv = [1.0,3.0]
    Ham.energies.NN = calc_E_NN( pw, strf, atoms.positions, atoms.Nspecies, atoms.atm2species, Zv)

    println("\nAfter calculating E_NN")
    println(Ham.energies)

    # states, need to be be set manually
    Nstates = 2
    Ham.focc = [2.0, 2.0]

    #
    KS_solve_Emin_PCG!( Ham, Nstates )
    println("\nAfter calling KS_solve_Emin_PCG:")
    println(Ham.energies)

end

@time test_main()
