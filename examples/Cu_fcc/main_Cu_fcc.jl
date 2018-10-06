function main()
    # Atoms
    atoms = Atoms(xyz_string_frac=
        """
        1

        Cu  0.0  0.0  0.0
        """, LatVecs = gen_lattice_fcc(3.61496*ANG2BOHR))

    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/Cu-q11.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5,
                       meshk=[8,8,8], extra_states=4 )
    println(Ham)

    #
    # Solve the KS problem
    #
    KS_solve_SCF!( Ham, betamix=0.2, mix_method="anderson", NiterMax=50, use_smearing=true )

end

