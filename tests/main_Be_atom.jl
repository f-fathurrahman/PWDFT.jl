using PWDFT

function test_main( ; method="SCF" )

    #
    # Atoms
    #
    atoms = init_atoms_xyz("Be.xyz")
    println(atoms)

    #
    # Initialize Hamiltonian
    #
    pspfiles = ["../pseudopotentials/pade_gth/Be-q2.gth"]
    LatVecs = 16.0*diagm( ones(3) )
    ecutwfc_Ry = 30.0
    Ham = PWHamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, LatVecs )

    println("sum V Ps loc = ", sum(Ham.potentials.Ps_loc))

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
        λ, v = KS_solve_Emin_PCG!( Ham, I_CG_BETA=4 )
        println("\nAfter calling KS_solve_Emin_PCG:")

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

@time test_main(method="Emin")

"""
Results from ABINIT: 30 Ry local pspot only
    Kinetic energy  =  1.50910017183361E+00
    Hartree energy  =  1.25207697685687E+00
    XC energy       = -7.86226407495056E-01
    Ewald energy    = -3.54662184935077E-01
    PspCore energy  = -1.39286404377276E-03
    Loc. psp. energy= -4.11236267761497E+00
    NL   psp  energy=  0.00000000000000E+00
    >>>>>>>>> Etotal= -2.49346698539839E+00

Results from ABINIT 30 Ry with NL pspot (the usual one)
    Kinetic energy  =  3.56342013836996E-01
    Hartree energy  =  3.50397710827049E-01
    XC energy       = -3.67782685733863E-01
    Ewald energy    = -3.54662184935077E-01
    PspCore energy  = -1.39286404377276E-03
    Loc. psp. energy= -1.21257971730573E+00
    NL   psp  energy=  2.38145519069942E-01
    >>>>>>>>> Etotal= -9.91532208284459E-01
"""
