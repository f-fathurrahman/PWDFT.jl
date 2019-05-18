function main( ; method="SCF" )
    # Atoms
    atoms = Atoms( xyz_string_frac=
        """
        4

        Zn      0.3333333   0.6666667   0.0000000
        Zn      0.6666667   0.3333333   0.5000000
        O       0.3333333   0.6666667   0.3450000
        O       0.6666667   0.3333333   0.8450000
        """, in_bohr=true,
        LatVecs = gen_lattice_hexagonal( 3.2495*ANG2BOHR, 5.2069*ANG2BOHR ) )

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Zn-q2.gth"),
                joinpath(DIR_PSP, "O-q6.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
    println(Ham)

    #
    # Solve the KS problem
    #
    if method == "SCF"
        KS_solve_SCF!( Ham, mix_method="rpulay" )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham, verbose=true )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    else
        error( @sprintf("ERROR: unknown method = %s", method) )
    end

end

