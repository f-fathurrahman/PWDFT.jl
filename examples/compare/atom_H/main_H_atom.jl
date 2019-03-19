function main( ; method="SCF" )
    # Atoms
    atoms = Atoms(xyz_string=
        """
        1

        H  0.0  0.0  0.0
        """, LatVecs=gen_lattice_sc(16.0))

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "H-q1.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc )

    println(Ham)

    #
    # Solve the KS problem
    #
    if method == "SCF"
        KS_solve_SCF!( Ham, mix_method="anderson", update_psi="LOBPCG" )

    elseif method == "CheFSI"
        KS_solve_SCF!( Ham, update_psi="CheFSI", betamix=0.5 )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    elseif method == "TRDCM"
        KS_solve_TRDCM!( Ham, NiterMax=50 )

    else
        error( @sprintf("Unknown method %s", method) )
    end

end
