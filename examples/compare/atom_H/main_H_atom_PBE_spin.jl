function main( ; method="SCF" )
    # Atoms
    atoms = Atoms( xyz_string=
        """
        1

        H  0.0  0.0  0.0
        """, LatVecs = gen_lattice_sc(16.0))

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PWDFT, "pseudopotentials", "pbe_gth", "H-q1.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, xcfunc="PBE",
                       Nspin=2, extra_states=0 )
    println(Ham)

    #
    # Solve the KS problem
    #
    if method == "SCF"
        KS_solve_SCF!( Ham, betamix=0.1, update_psi="LOBPCG", mix_method="anderson" )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham, verbose=true )

    elseif method == "DCM"
        KS_solve_DCM!( Ham )

    else
        error( @sprintf("ERROR: unknown method = %s", method) )
    end

end

