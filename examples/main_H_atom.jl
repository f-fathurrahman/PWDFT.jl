using PWDFT

function test_main( ; method="SCF" )

    #
    # Atoms
    #
    atoms = init_atoms_xyz("H.xyz")
    println(atoms)

    #
    # Initialize Hamiltonian
    #
    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth"]
    LatVecs = 16.0*diagm( ones(3) )
    ecutwfc_Ry = 30.0
    Ham = PWHamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, LatVecs )

    println("sum V Ps loc = ", sum(Ham.potentials.Ps_loc))

    #
    # calculate E_NN
    #
    Ham.energies.NN = calc_E_NN( Ham.pw, atoms, [1.0] )

    if method == "SCF"
        λ, v = KS_solve_SCF!( Ham )

    elseif method == "Emin"
        λ, v = KS_solve_Emin_PCG!( Ham, verbose=false )

    elseif method == "DCM"
        λ, v = KS_solve_DCM!( Ham )

    else
        println("ERROR: unknow method = ", method)
    end

    Nstates = Ham.electrons.Nstates
    println("\nEigenvalues")
    for ist = 1:Nstates
        @printf("%8d  %18.10f = %18.10f eV\n", ist, λ[ist], λ[ist]*Ry2eV*2)
    end
    println("\nTotal energy components")
    println(Ham.energies)

end

#@time test_main(method="Emin")
#@time test_main(method="Emin")
#@time test_main(method="DCM")
#@time test_main(method="SCF")
