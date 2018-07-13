using Printf
using PWDFT

function test_main( ; method="SCF" )

    # Atoms
    atoms = init_atoms_xyz("../../structures/Li2.xyz")
    atoms.LatVecs = gen_lattice_cubic(16.0)
    println(atoms)

    # Initialize Hamiltonian
    ecutwfc_Ry = 30.0
    pspfiles = ["../../pseudopotentials/pade_gth/Li-q1.gth"]
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
@time test_main(method="DCM")


#=
Result from ABINIT: (30 Ry):
    Kinetic energy  =  2.16029548905902E-01
    Hartree energy  =  1.78287277495300E-01
    XC energy       = -2.76638068416253E-01
    Ewald energy    = -2.19685854008068E-02
    PspCore energy  = -3.96586964629223E-03
    Loc. psp. energy= -6.14021970150543E-01
    NL   psp  energy=  1.48929934899351E-01
    >>>>>>>>> Etotal= -3.73347732313342E-01

=#
