using PWDFT

function test_main( ; method="SCF" )

    # Atoms
    atoms = init_atoms_xyz_string(
        """
        1

        C  0.0  0.0  0.0
        """
    )
    atoms.LatVecs = 16.0*diagm( ones(3) )
    println(atoms)

    # Initialize Hamiltonian
    ecutwfc_Ry = 40.0
    Ham = PWHamiltonian( atoms, ecutwfc_Ry*0.5 )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

    if method == "SCF"
        # Example of manually setting occupation numbers
        Ham.electrons.Nstates = 4
        Ham.electrons.Focc = [2.0, 4.0/3, 4.0/3, 4.0/3]
        KS_solve_SCF!( Ham )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    else
        println("ERROR: unknow method = ", method)

    end

    Nstates = Ham.electrons.Nstates
    ebands = Ham.electrons.ebands
    
    println("\nBand energies:")
    for ist = 1:Nstates
        @printf("%8d  %18.10f = %18.10f eV\n", ist, ebands[ist], ebands[ist]*Ry2eV*2)
    end
    
    println("\nTotal energy components")
    println(Ham.energies)

end

@time test_main(method="Emin")
@time test_main(method="SCF")
@time test_main(method="DCM")