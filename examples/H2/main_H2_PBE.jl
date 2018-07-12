using Printf
using PWDFT

function test_main( ; method="SCF" )

    # Atoms
    atoms = init_atoms_xyz("../../structures/H2.xyz")
    atoms.LatVecs = gen_lattice_sc(16.0)
    println(atoms)

    # Initialize Hamiltonian
    pspfiles = ["../../pseudopotentials/pbe_gth/H-q1.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, xcfunc="PBE" )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

    if method == "SCF"
        KS_solve_SCF!( Ham )

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

