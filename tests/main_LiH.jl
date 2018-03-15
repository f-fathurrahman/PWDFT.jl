using PWDFT

function test_main( ; method="SCF" )

    # Atoms
    atoms = init_atoms_xyz("LiH.xyz")
    println(atoms)

    # Initialize Hamiltonian
    LatVecs = 16.0*diagm( ones(3) )
    ecutwfc_Ry = 30.0
    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth",
                "../pseudopotentials/pade_gth/Li-q1.gth"]
    Ham = PWHamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, LatVecs )

    # calculate E_NN
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
        @printf("%8d  %18.10f\n", ist, λ[ist])
    end
    println("\nTotal energy components")
    println(Ham.energies)

end

@time test_main(method="Emin")
@time test_main(method="SCF")

"""
ABINIT result (30 Ry) Li-q1, H-q1
    Kinetic energy  =  6.12112937548124E-01
    Hartree energy  =  5.79304971112604E-01
    XC energy       = -4.78209829711666E-01
    Ewald energy    = -2.19685854008068E-02
    PspCore energy  = -1.98356853564844E-03
    Loc. psp. energy= -1.49966016707742E+00
    NL   psp  energy=  3.56955014824771E-02
    >>>>>>>>> Etotal= -7.74708740582338E-01

Total energy components
    Kinetic    energy:       0.6039873921
    Ps_loc     energy:      -1.4895615714
    Ps_nloc    energy:       0.0358362128
    Hartree    energy:       0.5766691587
    XC         energy:      -0.4770612342
    -------------------------------------
    Electronic energy:      -0.7501300421
    NN         energy:      -0.0219685854
    -------------------------------------
    Total      energy:      -0.7720986275
"""

