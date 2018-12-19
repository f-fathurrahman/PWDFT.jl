function main( ; method="SCF" )
    # Atoms
    atoms = Atoms( xyz_string_frac=
        """
        2

        Ga  0.0   0.0   0.0
        As  0.25  0.25  0.25
        """, in_bohr=true, LatVecs = gen_lattice_fcc(10.6839444516) )

    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/Ga-q3.gth",
                "../pseudopotentials/pade_gth/As-q5.gth"]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
    println(Ham)

    #
    # Solve the KS problem
    #
    if method == "SCF"
        KS_solve_SCF!( Ham, mix_method="anderson", betamix=0.2 )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham, verbose=true )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    else
        error( @sprintf("ERROR: unknown method = %s", method) )
    end

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
