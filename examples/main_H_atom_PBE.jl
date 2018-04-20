using PWDFT

function test_main( ; method="SCF" )

    #
    # Atoms
    #
    atoms = init_atoms_xyz("../structures/H.xyz")
    println(atoms)

    #
    # Initialize Hamiltonian
    #
    pspfiles = ["../pseudopotentials/pbe_gth/H-q1.gth"]
    LatVecs = 16.0*diagm( ones(3) )
    ecutwfc_Ry = 30.0
    Ham = PWHamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, LatVecs, xcfunc="PBE" )

    println("sum V Ps loc = ", sum(Ham.potentials.Ps_loc))

    #
    # calculate E_NN
    #
    Ham.energies.NN = calc_E_NN( Ham.pw, atoms, [1.0] )

    println("\nAfter calculating E_NN")
    println(Ham.energies)

    if method == "SCF"
        λ, v = KS_solve_SCF!( Ham )
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
        @printf("%8d  %18.10f = %18.10f eV\n", ist, λ[ist], λ[ist]*Ry2eV*2)
    end
    println("\nTotal energy components")
    println(Ham.energies)

end

@time test_main(method="Emin")
@time test_main(method="DCM")
@time test_main(method="SCF")

"""
ABINIT result (30 Ry)
    Kinetic energy  =  4.23404570165078E-01
    Hartree energy  =  1.98597587024920E-01
    XC energy       = -2.47350993458082E-01
    Ewald energy    = -8.86655462337693E-02
    PspCore energy  = -3.32190895798649E-07
    Loc. psp. energy= -7.38845265046175E-01
    NL   psp  energy=  0.00000000000000E+00
    >>>>>>>>> Etotal= -4.52859979738924E-01

Kinetic    energy:       0.3963474426
Ps_loc     energy:      -0.7091320512
Ps_nloc    energy:       0.0000000000
Hartree    energy:       0.1876606352
XC         energy:      -0.2384950342
-------------------------------------
Electronic energy:      -0.3636190076
NN         energy:      -0.0886655462
-------------------------------------
Total      energy:      -0.4522845539
"""
