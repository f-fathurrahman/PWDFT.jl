function main()
    # Atoms
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_fcc(3.9231*ANG2BOHR))

    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/Pt-q18.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, xcfunc="PBE",
                        meshk=[3,3,3], extra_states=4 )
    println(Ham)

    #
    # Solve the KS problem
    #
    KS_solve_SCF!( Ham, mix_method="anderson", betamix=0.2, use_smearing=true )

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

Kinetic    energy:      23.5372109366
Ps_loc     energy:      -5.3317439583
Ps_nloc    energy:     -13.8733354170
Hartree    energy:       3.2345191391
XC         energy:      -5.1246108125
-------------------------------------
Electronic energy:       2.4420398879
NN         energy:     -31.1867538034
-------------------------------------
Total      energy:     -28.7447139155
=#
