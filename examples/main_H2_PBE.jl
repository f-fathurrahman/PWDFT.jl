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
    pspfiles = ["../pseudopotentials/pbe_gth/H-q1.gth"]
    LatVecs = 16.0*diagm( ones(3) )
    ecutwfc_Ry = 30.0
    Ham = PWHamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, LatVecs, xcfunc="PBE" )


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
    Kinetic energy  =  1.03846627435554E+00
    Hartree energy  =  9.16452467490984E-01
    XC energy       = -6.66546679987766E-01
    Ewald energy    =  3.13170052325859E-01
    PspCore energy  = -1.32876358319460E-06
    Loc. psp. energy= -2.74851348731632E+00
    NL   psp  energy=  0.00000000000000E+00
    >>>>>>>>> Etotal= -1.14697270189529E+00

Kinetic    energy:       1.0022039947
Ps_loc     energy:      -2.7005240671
Ps_nloc    energy:       0.0000000000
Hartree    energy:       0.8935857782
XC         energy:      -0.6547008647
-------------------------------------
Electronic energy:      -1.4594351588
NN         energy:       0.3131700523
-------------------------------------
Total      energy:      -1.1462651065
"""

