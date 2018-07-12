using Printf
using PWDFT

function test_main( ; method="SCF" )

    # Atoms
    atoms = init_atoms_xyz("../../structures/N2.xyz")
    atoms.LatVecs = gen_lattice_cubic(16.0)
    println(atoms)

    # Initialize Hamiltonian
    ecutwfc_Ry = 30.0
    pspfiles = ["../../pseudopotentials/pade_gth/N-q5.gth"]
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5 )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

    if method == "SCF"
        KS_solve_SCF!( Ham, mix_method="anderson", β=0.2 )

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

@time test_main(method="Emin")
@time test_main(method="SCF")
@time test_main(method="DCM")


#=
Kinetic    energy:      11.5329676573
Ps_loc     energy:     -47.3445263895
Ps_nloc    energy:       2.3268197661
Hartree    energy:      17.1133647701
XC         energy:      -4.5072691732
-------------------------------------
Electronic energy:     -20.8786433691
NN         energy:       1.7906353998
-------------------------------------
Total      energy:     -19.0880079692

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
=#
