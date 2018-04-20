using PWDFT

function test_main( ; method="SCF" )

    #
    # Atoms
    #
    atoms = init_atoms_xyz("../structures/H2.xyz")
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
        println("\nAfter calling KS_solve_SCF:")

    elseif method == "ChebySCF"
        λ, v = KS_solve_SCF!( Ham, update_psi="CheFSI" )
        println("\nAfter calling KS_solve_SCF:")

    elseif method == "Emin"
        λ, v = KS_solve_Emin_PCG!( Ham )
        println("\nAfter calling KS_solve_Emin_PCG:")

    elseif method == "DCM"
        λ, v = KS_solve_DCM!( Ham )
        println("\nAfter calling KS_solve_Emin_DCM:")

    else
        println("ERROR: unknow method = ", method)
    end

    Nstates = Ham.electrons.Nstates
    println("\nEigenvalues")
    for ist = 1:Nstates
        @printf("%8d  %18.10f\n", ist, λ[ist])
    end
    println("\nTotal energy components")
    println(Ham.energies)

end

#@time test_main(method="ChebySCF")
@time test_main(method="Emin")
@time test_main(method="DCM")
#@time test_main(method="SCF")


"""
For 30 Ry

ABINIT result:
    Kinetic energy  =  1.00939508961180E+00
    Hartree energy  =  9.01100534759225E-01
    XC energy       = -6.30598700701720E-01
    Ewald energy    =  3.13170052325859E-01
    PspCore energy  = -1.26742500464741E-06
    Loc. psp. energy= -2.71195127559760E+00
    NL   psp  energy=  0.00000000000000E+00
    >>>>>>>>> Etotal= -1.11888556702743E+00


PWSCF result:
!    total energy              =      -2.23948549 Ry  = -1.119742745 Ha
     one-electron contribution =      -3.40513308 Ry
     hartree contribution      =       1.80222400 Ry
     xc contribution           =      -1.26291639 Ry
     ewald contribution        =       0.62633998 Ry  =  0.31316999 Ha
"""

