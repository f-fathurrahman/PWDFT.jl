function main( ; method="SCF" )
    # Atoms
    atoms = Atoms( xyz_string=
        """
        1

        H  0.0  0.0  0.0
        """, LatVecs = gen_lattice_sc(16.0))

    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, Nspin=2, extra_states=1 )
    println(Ham)

    #
    # Solve the KS problem
    #
    if method == "SCF"
        KS_solve_SCF!( Ham, mix_method="anderson" )

    elseif method == "CheFSI"
        KS_solve_SCF!( Ham, update_psi="CheFSI", betamix=0.5 )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham, verbose=true )

    elseif method == "DCM"
        KS_solve_DCM!( Ham )

    else
        error( @sprintf("ERROR: unknown method = %s", method) )
    end

end
