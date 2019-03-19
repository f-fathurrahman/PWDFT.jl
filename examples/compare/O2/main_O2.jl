function main( ; method="SCF" )

    # Atoms
    atoms = Atoms( xyz_file=joinpath(DIR_STRUCTURES, "O2.xyz"),
                   LatVecs = gen_lattice_cubic(16.0) )

    # Initialize Hamiltonian
    ecutwfc = 15.0
    pspfiles = [joinpath(DIR_PSP, "O-q6.gth")]
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc )
    println(Ham)

    if method == "SCF"
        KS_solve_SCF!( Ham, mix_method="anderson", betamix=0.2 )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    else
        error( @sprintf("Unknown method = %s", method) )
    end

end
