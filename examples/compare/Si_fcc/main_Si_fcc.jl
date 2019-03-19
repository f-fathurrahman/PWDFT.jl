function main( ; method="SCF" )
    # Atoms
    atoms = Atoms(xyz_string_frac=
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true, LatVecs=gen_lattice_fcc(10.2631))

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "Si-q4.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
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
    Kinetic energy  =  3.21045499991613E+00
    Hartree energy  =  5.76194835366404E-01
    XC energy       = -2.41016569027563E+00
    Ewald energy    = -8.39792740071415E+00
    PspCore energy  = -2.94625629171302E-01
    Loc. psp. energy= -2.17613632002650E+00
    NL   psp  energy=  1.58119195731965E+00
    >>>>>>>>> Etotal= -7.91101324758539E+00

PWDFT.jl:
Kinetic    energy:       3.2107141913
Ps_loc     energy:      -2.1756827897
Ps_nloc    energy:       1.5804698683
Hartree    energy:       0.5829626011
XC         energy:      -2.4156830816
PspCore    energy:      -0.2946256268
-------------------------------------
Electronic energy:       0.4881551626
NN         energy:      -8.3979274007
-------------------------------------
Total      energy:      -7.9097722381


QE:
one-elec   energy:       2.3207092600
Hartree    energy:       0.5764842900
XC         energy:      -2.4102748100
-------------------------------------
Electronic energy:       0.4869187400
NN         energy:      -8.3979274200
-------------------------------------
Total      energy:      -7.9110086800
=#
