push!(LOAD_PATH, "../../src")

using Printf
using PWDFT

function test_main( ; method="SCF" )

    # Atoms
    atoms = init_atoms_xyz("../../structures/He.xyz")
    atoms.LatVecs = gen_lattice_cubic(16.0)
    println(atoms)

    # Initialize Hamiltonian
    pspfiles = ["../../pseudopotentials/pade_gth/He-q2.gth"]
    ecutwfc_Ry = 30.0
    Ham = PWHamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5 )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

    if method == "SCF"
        KS_solve_SCF!( Ham )
        println("\nAfter calling KS_solve_SCF:")

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham )
        println("\nAfter calling KS_solve_Emin_PCG:")

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )
        println("\nAfter calling KS_solve_Emin_DCM:")

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

@time test_main(method="Emin")
@time test_main(method="DCM")
@time test_main(method="SCF")

"""
ABINIT (result):
    Kinetic energy  =  2.18802875961751E+00
    Hartree energy  =  1.50108932871126E+00
    XC energy       = -9.02904237573505E-01
    Ewald energy    = -3.54662184935077E-01
    PspCore energy  = -1.69161860064652E-06
    Loc. psp. energy= -5.08088416764117E+00
    NL   psp  energy=  0.00000000000000E+00
    >>>>>>>>> Etotal= -2.64933419343958E+00

Total energy components
Kinetic    energy:       2.1883489615
Ps_loc     energy:      -5.0813492117
Ps_nloc    energy:       0.0000000000
Hartree    energy:       1.5013376586
XC         energy:      -0.9038956763
-------------------------------------
Electronic energy:      -2.2955582679
NN         energy:      -0.3546621849
-------------------------------------
Total      energy:      -2.6502204529

"""
