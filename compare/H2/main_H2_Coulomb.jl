function main( ; method="SCF" )

    # Atoms
    atoms = Atoms( xyz_file="../structures/H2.xyz",
                   LatVecs=gen_lattice_sc(16.0) )

    # Initialize Hamiltonian
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, ecutwfc )
    println(Ham)

    if method == "SCF"
        KS_solve_SCF!( Ham )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    else
        error( @sprintf("ERROR: unknown method = %s", method) )
    end

end
