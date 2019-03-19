function main()
    # Atoms
    atoms = Atoms( xyz_string_frac=
        """
        1

        Cu  0.0  0.0  0.0
        """, LatVecs = gen_lattice_fcc(3.61496*ANG2BOHR) )

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Cu-q11.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, Nspin=2, xcfunc="VWN",
                         meshk=[3,3,3], extra_states=4 )

    #
    # Solve the KS problem
    #
    KS_solve_SCF!( Ham, betamix=0.5, mix_method="anderson", use_smearing=true )

end
