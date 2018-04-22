using PWDFT

function test_main( ; method="SCF" )

    # Atoms
    atoms = init_atoms_xyz("../structures/CO.xyz")
    println(atoms)

    # Initialize Hamiltonian
    LatVecs = 16.0*diagm( ones(3) )
    ecutwfc_Ry = 30.0
    pspfiles = ["../pseudopotentials/pade_gth/C-q4.gth",
                "../pseudopotentials/pade_gth/O-q6.gth"]
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

@time test_main(method="DCM") # diverges
#@time test_main(method="Emin")
#@time test_main(method="SCF")


"""
    Kinetic energy  =  1.27114662836776E+01
    Hartree energy  =  1.87542476647748E+01
    XC energy       = -4.73571914335774E+00
    Ewald energy    =  2.26310215311565E+00
    PspCore energy  = -2.54693990241102E-04
    Loc. psp. energy= -5.21197227228070E+01
    NL   psp  energy=  2.59183669792393E+00
    >>>>>>>>> Etotal= -2.05350437606631E+01

Kinetic    energy:      12.7122820382
Ps_loc     energy:     -52.1220967254
Ps_nloc    energy:       2.5918770838
Hartree    energy:      18.7558153675
XC         energy:      -4.7410160310
-------------------------------------
Electronic energy:     -22.8031382669
NN         energy:       2.2631029868
-------------------------------------
Total      energy:     -20.5400352801
"""
