function main()
    # Atoms
    atoms = Atoms( xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true,
        LatVecs = gen_lattice_fcc(10.2631) )

    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/Si-q4.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, xcfunc="VWN",
                       meshk=[3,3,3], Nspin=2, extra_states=0 )
    println(Ham)

    #
    # Solve the KS problem
    #
    KS_solve_SCF!( Ham, mix_method="anderson" )

end

#=
No smearing:

Total energy components

Kinetic    energy:       3.1744863175
Ps_loc     energy:      -2.3423254958
Ps_nloc    energy:       1.5694366997
Hartree    energy:       0.5730118848
XC         energy:      -2.4296659652
-------------------------------------
Electronic energy:       0.5449434409
NN         energy:      -8.3979274007
-------------------------------------
Total      energy:      -7.8529839598

--------------

With smearing 0.01 Ha
Total energy components

Kinetic    energy:       3.2023421562
Ps_loc     energy:      -2.2762398006
Ps_nloc    energy:       1.5915987675
Hartree    energy:       0.5322060942
XC         energy:      -2.2975477892
-------------------------------------
Electronic energy:       0.7523594282
NN         energy:      -8.3979274007
-------------------------------------
Total      energy:      -7.6455679725

------------

ABINIT

Kinetic energy  =  3.18425700948638E+00
Hartree energy  =  5.65422931914669E-01
XC energy       = -2.42294101779165E+00
Ewald energy    = -8.39792740071415E+00
PspCore energy  = -2.09890966765802E-01
Loc. psp. energy= -2.12992717173527E+00
NL   psp  energy=  1.55983881823907E+00
>>>>> Internal E= -7.85116779736675E+00

-kT*entropy     = -4.52029451264488E-03
>>>>>>>>> Etotal= -7.85568809187940E+00

PWSCF

!    total energy              =     -15.70835819 Ry = -7.854179095
     Harris-Foulkes estimate   =     -15.70835757 Ry
     estimated scf accuracy    <       0.00000013 Ry
     one-electron contribution =       4.79575036 Ry
     hartree contribution      =       1.14632766 Ry
     xc contribution           =      -4.85295131 Ry
     ewald contribution        =     -16.79585054 Ry
     smearing contrib. (-TS)   =      -0.00163437 Ry

=#
