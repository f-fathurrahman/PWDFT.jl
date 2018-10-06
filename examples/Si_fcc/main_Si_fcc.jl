function main( ; method="SCF" )
    # Atoms
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, LatVecs=gen_lattice_fcc(10.2631))

    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/Si-q4.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, meshk=[3,3,3] )
    println(Ham)

    #
    # Solve the KS problem
    #
    if method == "SCF"
        KS_solve_SCF!( Ham, mix_method="anderson" )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham, verbose=true )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    elseif method == "TRDCM"
        KS_solve_TRDCM!( Ham, NiterMax=50 )

    else
        error( @sprintf("Unknown method %s", method) )
    end
    
end

#=
FFT grid = (27,27,27)
!    total energy              =     -15.82591729 Ry = -7.912958645
     one-electron contribution =       4.64239409 Ry
     hartree contribution      =       1.15137054 Ry
     xc contribution           =      -4.82383139 Ry
     ewald contribution        =     -16.79585054 Ry = -8.39792527

Kinetic energy  =  3.20979169051178E+00
Hartree energy  =  5.75704916874837E-01
XC energy       = -2.40756461272468E+00
Ewald energy    = -8.39792740071415E+00
PspCore energy  = -2.94625629171302E-01
Loc. psp. energy= -2.17494660308169E+00
NL   psp  energy=  1.58096345360260E+00
>>>>>>>>> Etotal= -7.90860418470260E+00

Kinetic    energy:       3.2108017168
Ps_loc     energy:      -2.1760104813
Ps_nloc    energy:       1.5806486758
Hartree    energy:       0.5830565701
XC         energy:      -2.4157170879
-------------------------------------
Electronic energy:       0.7827793934
NN         energy:      -8.3979258900
-------------------------------------
Total      energy:      -7.6151464966
PspCore ene = -0.29462546778533705

TotEne + PspCore = -7.909771964367133


FFT grid = (27,27,27)
Kinetic    energy:       3.2107925351
Ps_loc     energy:      -2.4706316946
Ps_nloc    energy:       1.5806564220
Hartree    energy:       0.5830491953
XC         energy:      -2.4157133347
-------------------------------------
Electronic energy:       0.4881531232
NN         energy:      -8.3979258900
-------------------------------------
Total      energy:      -7.9097727668
=#
