using PWDFT

function test_main( ; method="SCF" )
    const LatVecs = 16.0*diagm( ones(3) )
    ecutwfc_Ry = 40.0
    pw = PWGrid( ecutwfc_Ry*0.5, LatVecs )
    println(pw)

    #
    # Atoms
    #
    atoms = init_atoms_xyz("LiH.xyz")
    println(atoms)

    #
    # Structure factor
    #
    strf = calc_strfact( atoms, pw )

    #
    # Initialize Hamiltonian
    #
    Ham = PWHamiltonian(pw)
    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth",
                "../pseudopotentials/pade_gth/Li-q3.gth"]
    update!(Ham, atoms, strf, pspfiles)
    println("sum V Ps loc = ", sum(Ham.potentials.Ps_loc))

    #
    # calculate E_NN
    #
    Zv = [1.0,3.0]
    Ham.energies.NN = calc_E_NN( pw, strf, atoms.positions, atoms.Nspecies, atoms.atm2species, Zv)

    println("\nAfter calculating E_NN")
    println(Ham.energies)

    # states, need to be be set manually
    Nstates = 2
    Ham.focc = [2.0, 2.0]

    if method == "SCF"
        λ, v = KS_solve_SCF!( Ham, Nstates, update_psi="CheFSI" )
        println("\nAfter calling KS_solve_SCF:")

    elseif method == "Emin"
        λ, v = KS_solve_Emin_PCG!( Ham, Nstates, I_CG_BETA=4 )
        println("\nAfter calling KS_solve_Emin_PCG:")

    else
        println("ERROR: unknow method = ", method)
    end

    println("\nEigenvalues")
    for ist = 1:Nstates
        @printf("%8d  %18.10f\n", ist, λ[ist])
    end
    println("\nTotal energy components")
    println(Ham.energies)

end

@time test_main(method="Emin")
#@time test_main(method="SCF")
