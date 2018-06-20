using Printf
using PWDFT

function test_main( newpos ; method="SCF" )
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        1

        H  0.0  0.0  0.0
        """)
    atoms.LatVecs = gen_lattice_sc(16.0)
    atoms.positions[:,1] = newpos  # manually set the position
    println(atoms)

    # Initialize Hamiltonian
    pspfiles = ["../../pseudopotentials/pade_gth/H-q1.gth"]
    ecutwfc_Ry = 30.0
    Ham = PWHamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5 )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

    #
    # Solve the KS problem
    #
    if method == "SCF"
        KS_solve_SCF!( Ham, mix_method="anderson", ETOT_CONV_THR=1e-6 )

    elseif method == "CheFSI"
        KS_solve_SCF!( Ham, update_psi="CheFSI", β=0.5 )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham, verbose=true )

    elseif method == "DCM"
        KS_solve_DCM!( Ham )

    else
        println("ERROR: unknown method = ", method)
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

@time test_main([0.0, 0.0, 0.0], method="Emin")
@time test_main([0.0, 0.0, 0.05], method="Emin")
@time test_main([0.0, 0.0, 0.1], method="Emin")
@time test_main([0.0, 0.0, 0.15], method="Emin")
@time test_main([8.0, 8.0, 8.0], method="Emin")