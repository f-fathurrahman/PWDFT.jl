using PWDFT

function test_main( ; method="SCF" )

    #
    # Atoms
    #
    atoms = init_atoms_xyz("Li2.xyz")
    println(atoms)

    #
    # Initialize Hamiltonian
    #
    LatVecs = 16.0*diagm( ones(3) )
    ecutwfc_Ry = 30.0
    pspfiles = ["../pseudopotentials/pade_gth/Li-q1.gth"]
    Ham = PWHamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, LatVecs )

    #
    # calculate E_NN
    #
    Zvals = get_Zvals( Ham.pspots )
    Ham.energies.NN = calc_E_NN( Ham.pw, atoms, Zvals )

    println("\nAfter calculating E_NN")
    println(Ham.energies)

    if method == "SCF"
        λ, v = KS_solve_SCF!( Ham )
        println("\nAfter calling KS_solve_SCF:")

    elseif method == "Emin"
        λ, v = KS_solve_Emin_PCG!( Ham )
        println("\nAfter calling KS_solve_Emin_PCG:")

    else
        println("ERROR: unknow method = ", method)
    end

    Nstates = Ham.electrons.Nstates
    println("\nEigenvalues")
    for ist = 1:Nstates
        @printf("%8d  %18.10f = %18.10f eV\n", ist, λ[ist], λ[ist]*2.0*Ry2eV)
    end
    println("\nTotal energy components")
    println(Ham.energies)

end

@time test_main(method="Emin")
#@time test_main(method="SCF")


"""
Result from ABINIT: (30 Ry):
    Kinetic energy  =  2.16029548905902E-01
    Hartree energy  =  1.78287277495300E-01
    XC energy       = -2.76638068416253E-01
    Ewald energy    = -2.19685854008068E-02
    PspCore energy  = -3.96586964629223E-03
    Loc. psp. energy= -6.14021970150543E-01
    NL   psp  energy=  1.48929934899351E-01
    >>>>>>>>> Etotal= -3.73347732313342E-01

"""
