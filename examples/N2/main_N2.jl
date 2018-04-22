using PWDFT

function test_main( ; method="SCF" )

    #
    # Atoms
    #
    atoms = init_atoms_xyz("../structures/N2.xyz")
    println(atoms)

    #
    # Initialize Hamiltonian
    #
    LatVecs = 16.0*diagm( ones(3) )
    ecutwfc_Ry = 30.0
    pspfiles = ["../pseudopotentials/pade_gth/N-q5.gth"]
    Ham = PWHamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, LatVecs )

    @printf("\nsum V Ps loc = %18.10f\n", sum(Ham.potentials.Ps_loc))

    #
    # calculate E_NN
    #
    Zvals = get_Zvals( Ham.pspots )
    Ham.energies.NN = calc_E_NN( Ham.pw, atoms, Zvals )

    if method == "SCF"
        λ, v = KS_solve_SCF!( Ham, β=0.2 )
        println("\nAfter calling KS_solve_SCF:")
    
    elseif method == "Emin"
        λ, v = KS_solve_Emin_PCG!( Ham )
        println("\nAfter calling KS_solve_Emin_PCG:")
 
    elseif method == "DCM"
        λ, v = KS_solve_DCM!( Ham )
        println("\nAfter calling KS_solve_DCM:")

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
@time test_main(method="DCM")
#@time test_main(method="SCF")


"""
ABINIT result for 30 Ry:
    Kinetic energy  =  1.15309586059089E+01
    Hartree energy  =  1.71105433424343E+01
    XC energy       = -4.50162677653213E+00
    Ewald energy    =  1.79063539962452E+00
    PspCore energy  = -7.02139897582120E-05
    Loc. psp. energy= -4.73400706683002E+01
    NL   psp  energy=  2.32674629031602E+00
    >>>>>>>>> Etotal= -1.90828840205384E+01


PWSCF result for 30 Ry
!    total energy              =     -38.17145300 Ry = -19.0857265 Ha
     one-electron contribution =     -66.96020106 Ry
     hartree contribution      =      34.21524310 Ry
     xc contribution           =      -9.00776423 Ry
     ewald contribution        =       3.58126919 Ry =   1.790634595 Ha
"""
