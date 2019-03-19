function main()
    # Atoms
    atoms = Atoms( xyz_string_frac=
        """
        1

        Ni  0.0  0.0  0.0
        """, in_bohr=true,
        LatVecs = gen_lattice_fcc(6.65914911201) )

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Ni-q18.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, xcfunc="LDA",
                       Nspin=2, meshk=[3,3,3], extra_states=4 )
    println(Ham)

    #
    # Solve the KS problem
    #
    KS_solve_SCF!( Ham, mix_method="rpulay", betamix=0.2, use_smearing=true, kT=0.01 )
    
end

