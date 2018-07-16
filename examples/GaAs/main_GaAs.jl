function main( ; method="SCF" )
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        2

        Ga  0.0   0.0   0.0
        As  0.25  0.25  0.25
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(5.6537*ANG2BOHR)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)
    write_xsf( "TEMP_GaAs.xsf", atoms )

    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/Ga-q3.gth",
                "../pseudopotentials/pade_gth/As-q5.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, meshk=[4,4,4], verbose=true )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

    #
    # Solve the KS problem
    #
    if method == "SCF"
        # FIXME: need a more effective ways to deal with this
        Ham.electrons = Electrons( atoms, Ham.pspots, Nstates=5,
                                   Nkpt=Ham.pw.gvecw.kpoints.Nkpt, Nstates_empty=1 )
        println(Ham.electrons)
        KS_solve_SCF!( Ham, mix_method="anderson", β=0.2 )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham, verbose=true )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    else
        println("ERROR: unknown method = ", method)
    end

    println("\nTotal energy components")
    println(Ham.energies)

    pspcore_ene = calc_PspCore_ene(atoms, Ham.pspots, Ham.pw.Ω)
    
    println("\nPspCore ene = ", pspcore_ene)
    println("\nTotEne + PspCore = ", pspcore_ene + Ham.energies.Total)


end

#=
!    total energy              =     -17.28495473 Ry = -8.642477365 Ha
     one-electron contribution =       2.73205378 Ry
     hartree contribution      =       1.63341388 Ry
     xc contribution           =      -4.80801106 Ry
     ewald contribution        =     -16.84241132 Ry = -4.21060283 Ha

Kinetic energy  =  3.26683000819830E+00
Hartree energy  =  8.17221804456118E-01
XC energy       = -2.40229498371500E+00
Ewald energy    = -8.42120578386954E+00
PspCore energy  =  3.78008143191317E-01
Loc. psp. energy= -3.13731037111346E+00
NL   psp  energy=  8.58159437985374E-01
>>>>>>>>> Etotal= -8.64059174486689E+00 

Kinetic    energy:       3.2350078422
Ps_loc     energy:      -2.7329323313
Ps_nloc    energy:       0.8610819832
Hartree    energy:       0.8024651116
XC         energy:      -2.3993127351
-------------------------------------
Electronic energy:      -0.2336901294
NN         energy:      -8.4212062785
-------------------------------------
Total      energy:      -8.6548964079
=#
