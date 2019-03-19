function main()
    # Atoms
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pd  0.0  0.0  0.0
        """, LatVecs=gen_lattice_fcc(7.35065658378))

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Pd-q10.gth")]
    ecutwfc = 40.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc,
                       meshk=[3,3,3], extra_states=4 )
    println(Ham)

    #
    # Solve the KS problem
    #
    KS_solve_SCF!( Ham, use_smearing=true, mix_method="rpulay", kT=0.01 )

end

#=
!    total energy              =     -57.24503621 Ry = -28.622518105 Ha
     one-electron contribution =       9.04437856 Ry
     hartree contribution      =       6.25916206 Ry
     xc contribution           =      -9.63113297 Ry
     ewald contribution        =     -62.37350289 Ry = -31.186751445 Ha
     smearing contrib. (-TS)   =      -0.54394096 Ry =  -0.27197048 Ha

    Kinetic energy  =  2.43064410898704E+01
    Hartree energy  =  3.12650553185174E+00
    XC energy       = -4.81248242615583E+00
    Ewald energy    = -3.11867519714552E+01
    PspCore energy  =  3.99722461197524E+00
    Loc. psp. energy= -9.40369691424283E+00
    NL   psp  energy= -1.44442917976521E+01
    >>>>> Internal E= -2.84170518758087E+01

    -kT*entropy     = -2.26701305153683E-02
    >>>>>>>>> Etotal= -2.84397220063240E+01

Kinetic    energy:      25.5438065335
Ps_loc     energy:      -9.8715679296
Ps_nloc    energy:     -16.0142752833
Hartree    energy:       3.4706813329
XC         energy:      -4.8667156881
PspCore    energy:       3.9972245314
-TS              :      -0.0222154525
-------------------------------------
Electronic energy:       2.2369380442
NN         energy:     -31.1867519715
-------------------------------------
Total free energy:     -28.9498139273
=#
