using PWDFT

function test_main( ; method="SCF" )

    #
    # Atoms
    #
    atoms = init_atoms_xyz("../structures/NH3.xyz")
    println(atoms)

    #
    # Initialize Hamiltonian
    #
    LatVecs = 16.0*diagm( ones(3) )
    ecutwfc_Ry = 30.0
    pspfiles = ["../pseudopotentials/pade_gth/N-q5.gth",
                "../pseudopotentials/pade_gth/H-q1.gth"]
    Ham = PWHamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, LatVecs )

    # calculate E_NN
    Zvals = get_Zvals( Ham.pspots )
    Ham.energies.NN = calc_E_NN( Ham.pw, atoms, Zvals )

    if method == "SCF"
        λ, v = KS_solve_SCF!( Ham, β=0.2 )
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
        @printf("%8d  %18.10f = %18.10f eV\n", ist, λ[ist], λ[ist]*Ry2eV*2)
    end
    println("\nTotal energy components")
    println(Ham.energies)

end

@time test_main(method="Emin")
@time test_main(method="SCF")

"""
ABINIT 30 Ry
Kinetic energy  =  8.06714954093561E+00
Hartree energy  =  1.18361035301009E+01
XC energy       = -3.43865606332895E+00
Ewald energy    =  3.13315248577189E+00
PspCore energy  = -3.56901459311693E-05
Loc. psp. energy= -3.18586718047665E+01
NL   psp  energy=  9.94081601344134E-01
>>>>>>>>> Etotal= -1.12668764000888E+01


Total energy components
Kinetic    energy:       8.0683940906
Ps_loc     energy:     -31.8614308498
Ps_nloc    energy:       0.9939896839
Hartree    energy:      11.8380732002
XC         energy:      -3.4429522512
-------------------------------------
Electronic energy:     -14.4039261263
NN         energy:       3.1331531467
-------------------------------------
Total      energy:     -11.2707729796
"""
