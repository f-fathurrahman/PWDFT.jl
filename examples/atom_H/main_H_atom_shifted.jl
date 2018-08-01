function main(; method="SCF")
    do_pos([0.0, 0.0, 0.0], "Emin")
    do_pos([0.0, 0.0, 0.05], "Emin")
    do_pos([0.0, 0.0, 0.1], "Emin")
    do_pos([0.0, 0.0, 0.15], "Emin")
    do_pos([8.0, 8.0, 8.0], "Emin")
end

function do_pos( newpos::Array{Float64,1}, method::String )
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
    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5 )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

    #
    # Solve the KS problem
    #
    if method == "SCF"
        KS_solve_SCF!( Ham, mix_method="anderson", ETOT_CONV_THR=1e-6 )

    elseif method == "CheFSI"
        KS_solve_SCF!( Ham, update_psi="CheFSI", betamix=0.5 )

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