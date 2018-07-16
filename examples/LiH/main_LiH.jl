function main( ; method="SCF" )

    # Atoms
    atoms = init_atoms_xyz("../structures/LiH.xyz")
    atoms.LatVecs = gen_lattice_cubic(16.0)
    println(atoms)

    # Initialize Hamiltonian
    ecutwfc_Ry = 30.0
    pspfiles = ["../pseudopotentials/pade_gth/Li-q3.gth",
                "../pseudopotentials/pade_gth/H-q1.gth"]
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5 )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

    if method == "SCF"
        KS_solve_SCF!( Ham, mix_method="anderson" )

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

#=
ABINIT result (30 Ry) Li-q1, H-q1
    Kinetic energy  =  6.03326166300386E-01
    Hartree energy  =  5.76057114250496E-01
    XC energy       = -4.76473283615080E-01
    Ewald energy    = -2.19685854008068E-02
    PspCore energy  = -1.98356853564844E-03
    Loc. psp. energy= -1.48661401427370E+00
    NL   psp  energy=  3.58557613147743E-02
    >>>>>>>>> Etotal= -7.71800409959577E-01

Total energy components (30 Ry) Li-q1, H-q1
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

-------------------------------------------------------

ABINIT (30 Ry), Li-q3, H-q1
    Kinetic energy  =  5.21627540905992E+00
    Hartree energy  =  3.60118849544695E+00
    XC energy       = -1.79213684515667E+00
    Ewald energy    = -4.20567941137495E-01
    PspCore energy  = -2.17841436121264E-05
    Loc. psp. energy= -1.36038024841280E+01
    NL   psp  energy=  0.00000000000000E+00
    >>>>>>>>> Etotal= -6.99906515005892E+00

Total energy components

Kinetic    energy:       5.2166027677
Ps_loc     energy:     -13.6047265039
Ps_nloc    energy:       0.0000000000
Hartree    energy:       3.6020560721
XC         energy:      -1.7936808540
-------------------------------------
Electronic energy:      -6.5797485182
NN         energy:      -0.4205679411
-------------------------------------
Total      energy:      -7.0003164593
----------------------------------------------------------

=#

