function main( ; method="SCF" )

    # Atoms
    atoms = Atoms( xyz_file="../structures/H2.xyz",
                   LatVecs = gen_lattice_sc(16.0) )

    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth"]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc )
    println(Ham)

    if method == "SCF"
        KS_solve_SCF!( Ham, mix_method="pulay" )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    elseif method == "TRDCM"
        KS_solve_TRDCM!( Ham, NiterMax=15 )

    else
        error( @sprintf("ERROR: unknown method = %s", method) )
    end

end

#=
For 30 Ry

ABINIT result:
    Kinetic energy  =  1.00939508961180E+00
    Hartree energy  =  9.01100534759225E-01
    XC energy       = -6.30598700701720E-01
    Ewald energy    =  3.13170052325859E-01
    PspCore energy  = -1.26742500464741E-06
    Loc. psp. energy= -2.71195127559760E+00
    NL   psp  energy=  0.00000000000000E+00
    >>>>>>>>> Etotal= -1.11888556702743E+00


PWSCF result:
!    total energy              =      -2.23948549 Ry  = -1.119742745 Ha
     one-electron contribution =      -3.40513308 Ry
     hartree contribution      =       1.80222400 Ry
     xc contribution           =      -1.26291639 Ry
     ewald contribution        =       0.62633998 Ry  =  0.31316999 Ha
=#

