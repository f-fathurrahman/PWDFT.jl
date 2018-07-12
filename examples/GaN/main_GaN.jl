using Printf
using PWDFT

function test_main( ; method="SCF" )
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        4

        Ga   0.333333333333333   0.666666666666667   0.000000000000000 
        Ga   0.666666666666667   0.333333333333333   0.500000000000000 
         N   0.333333333333333   0.666666666666667   0.385000000000000 
         N   0.666666666666667   0.333333333333333   0.885000000000000
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_hexagonal( 3.18*ANG2BOHR, 5.166*ANG2BOHR )
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)
    write_xsf( "TEMP_GaN.xsf", atoms )

    # Initialize Hamiltonian
    pspfiles = ["../../pseudopotentials/pade_gth/Ga-q3.gth",
                "../../pseudopotentials/pade_gth/N-q5.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, meshk=[3,3,3], verbose=true )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

    #
    # Solve the KS problem
    #
    if method == "SCF"
        KS_solve_SCF!( Ham )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham, verbose=true )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    else
        println("ERROR: unknow method = ", method)
    end

    println("\nTotal energy components")
    println(Ham.energies)

end

@time test_main(method="Emin")
@time test_main(method="SCF")
@time test_main(method="DCM")


#=
!    total energy              =     -47.51053311 Ry = -23.755266555 Ha
     one-electron contribution =      -1.86091588 Ry
     hartree contribution      =       9.80625212 Ry
     xc contribution           =     -13.10323777 Ry
     ewald contribution        =     -42.35263158 Ry = -21.17631579 Ha

    Kinetic energy  =  1.44378564168680E+01
    Hartree energy  =  4.90773882272804E+00
    XC energy       = -6.55396420607758E+00
    Ewald energy    = -2.11763161375599E+01
    PspCore energy  =  6.18064119181681E-01
    Loc. psp. energy= -1.83151350511424E+01
    NL   psp  energy=  2.32530222716897E+00
    >>>>>>>>> Etotal= -2.37564538088332E+01


Kinetic    energy:      14.4370546419
Ps_loc     energy:     -17.6965369131
Ps_nloc    energy:       2.3247884467
Hartree    energy:       4.9157426294
XC         energy:      -6.5597802583
-------------------------------------
Electronic energy:      -2.5787314534
NN         energy:     -21.1763173814
-------------------------------------
Total      energy:     -23.7550488349
=#
