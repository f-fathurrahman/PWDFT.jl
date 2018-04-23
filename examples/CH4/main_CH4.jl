using PWDFT

function test_main( ; method="SCF" )

    # Atoms
    atoms = init_atoms_xyz("../../structures/CH4.xyz")
    atoms.LatVecs = 16.0*diagm( ones(3) )    
    println(atoms)

    # Initialize Hamiltonian
    ecutwfc_Ry = 30.0
    pspfiles = ["../../pseudopotentials/pade_gth/C-q4.gth",
                "../../pseudopotentials/pade_gth/H-q1.gth"]
    Ham = PWHamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5 )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

    if method == "SCF"
        KS_solve_SCF!( Ham, β=0.2 )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    else
        println("ERROR: unknow method = ", method)
    end

    Nstates = Ham.electrons.Nstates
    ebands = Ham.electrons.ebands
    
    println("\nBand energies:")
    for ist = 1:Nstates
        @printf("%8d  %18.10f = %18.10f eV\n", ist, ebands[ist], ebands[ist]*Ry2eV*2)
    end
    
    println("\nTotal energy components")
    println(Ham.energies)

end

@time test_main(method="DCM") # diverges
@time test_main(method="Emin")
@time test_main(method="SCF")


"""
Kinetic energy  =  6.21988590075257E+00
Hartree energy  =  9.51102621859118E+00
XC energy       = -3.06030478665958E+00
Ewald energy    =  3.92782713134035E+00
PspCore energy  = -3.41579716509506E-04
Loc. psp. energy= -2.49780113532849E+01
NL   psp  energy=  5.00640903711708E-01
>>>>>>>>> Etotal= -7.87927756526521E+00

Kinetic    energy:       6.2219148018
Ps_loc     energy:     -24.9826270157
Ps_nloc    energy:       0.5007439033
Hartree    energy:       9.5137557131
XC         energy:      -3.0644609033
-------------------------------------
Electronic energy:     -11.8106735008
NN         energy:       3.9278278464
-------------------------------------
Total      energy:      -7.8828456544
"""
