using Printf
using PWDFT

function test_main( ; method="SCF" )

    # Atoms
    atoms = init_atoms_xyz("../../structures/CO2.xyz")
    atoms.LatVecs = gen_lattice_cubic(16.0)
    println(atoms)

    # Initialize Hamiltonian
    ecutwfc_Ry = 30.0
    pspfiles = ["../../pseudopotentials/pade_gth/C-q4.gth",
                "../../pseudopotentials/pade_gth/O-q6.gth"]
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5 )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

    if method == "SCF"
        KS_solve_SCF!( Ham )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    else
        println("ERROR: unknown method = ", method)
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

@time test_main(method="Emin")
@time test_main(method="SCF")
@time test_main(method="DCM") # diverges

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
