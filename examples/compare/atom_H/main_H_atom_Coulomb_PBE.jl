function main( ; method="SCF" )
    # Atoms
    atoms = Atoms(xyz_string=
        """
        1

        H  0.0  0.0  0.0
        """, LatVecs = gen_lattice_cubic(16.0))

    # Initialize Hamiltonian
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, ecutwfc, xcfunc="PBE" )
    println(Ham)

    #
    # Solve the KS problem
    #
    if method == "SCF"
        KS_solve_SCF!( Ham )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    elseif method == "TRDCM"
        KS_solve_TRDCM!( Ham )

    else
        error( @sprintf("Unknown method = %s", method) )
    end

end
