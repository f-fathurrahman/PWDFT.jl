function main()
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        1

        Ni  0.0  0.0  0.0
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(6.65914911201)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)

    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/Ni-q18.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, xcfunc="LDA",
                       Nspin=2, meshk=[3,3,3], verbose=true, extra_states=4 )
    #Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    #Ham.electrons.Focc[:,1:Nkpt] .= [1.0, 1.0, 1.0, 1.0, 1.0,
    #                                 1.0, 1.0, 1.0, 1.0, 1.0,
    #                                 0.0, 0.0, 0.0]
    #Ham.electrons.Focc[:,Nkpt+1:2*Nkpt] .= [1.0, 1.0, 1.0, 1.0, 1.0,
    #                                        1.0, 1.0, 1.0, 0.0, 0.0,
    #                                        0.0, 0.0, 0.0]



    #
    # Solve the KS problem
    #
    KS_solve_SCF!( Ham, mix_method="anderson", betamix=0.2, use_smearing=true )
    
    println("\nTotal energy components")
    println(Ham.energies)

end

