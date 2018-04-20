using PWDFT

function test_main( ; method="SCF" )

    # Atoms
    atoms = init_atoms_xyz("../structures/CO2.xyz")
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

#@time test_main(method="DCM") # diverges
@time test_main(method="Emin")
@time test_main(method="SCF")


"""
    Kinetic energy  =  2.23593227243830E+01
    Hartree energy  =  3.63151213069911E+01
    XC energy       = -8.05250269608862E+00
    Ewald energy    =  7.31427913787265E+00
    PspCore energy  = -1.52140135826873E-04
    Loc. psp. energy= -9.78903195421267E+01
    NL   psp  energy=  4.34788251847085E+00
    >>>>>>>>> Etotal= -3.56063686906335E+01

Kinetic    energy:      22.3603997337
Ps_loc     energy:     -97.8942867088
Ps_nloc    energy:       4.3479256820
Hartree    energy:      36.3182773651
XC         energy:      -8.0614719027
-------------------------------------
Electronic energy:     -42.9291558307
NN         energy:       7.3142812924
-------------------------------------
Total      energy:     -35.6148745383
"""
