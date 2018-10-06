function main( ; method="SCF" )

    # Atoms
    atoms = Atoms(xyz_file="../structures/NH3.xyz",
                  LatVecs=gen_lattice_cubic(16.0))

    # Initialize Hamiltonian
    ecutwfc_Ry = 30.0
    pspfiles = ["../pseudopotentials/pade_gth/N-q5.gth",
                "../pseudopotentials/pade_gth/H-q1.gth"]
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5 )
    println(Ham)

    if method == "SCF"
        KS_solve_SCF!( Ham )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    elseif method == "TRDCM"
        KS_solve_TRDCM!( Ham, NiterMax=50 )

    else
        error( @sprintf("Unknown method %s", method) )
    end

end

#=
ABINIT 30 Ry
    Kinetic energy  =  8.06862489714732E+00
    Hartree energy  =  1.18383732395827E+01
    XC energy       = -3.44300449580879E+00
    Ewald energy    =  3.13315248577189E+00
    PspCore energy  = -3.56901459311693E-05
    Loc. psp. energy= -3.18618998793985E+01
    NL   psp  energy=  9.94019129116458E-01
    >>>>>>>>> Etotal= -1.12707703137349E+01

Total energy components
Kinetic    energy:       8.0683940906
Ps_loc     energy:     -31.8614308498
Ps_nloc    energy:       0.9939896839
Hartree    energy:      11.8380732002
XC         energy:      -3.4429522512
-------------------------------------
Electronic energy:     -14.4039261263
NN         energy:       3.1331531467
-------------------------------------
Total      energy:     -11.2707729796
=#
