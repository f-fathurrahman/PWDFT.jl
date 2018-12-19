function main( ; method="SCF" )

    # Atoms
    atoms = Atoms( xyz_file="../structures/Li2.xyz",
                   LatVecs=gen_lattice_sc(16.0) )

    # Initialize Hamiltonian
    ecutwfc = 15.0
    pspfiles = ["../pseudopotentials/pade_gth/Li-q1.gth"]
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc )

    if method == "SCF"
        KS_solve_SCF!( Ham )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    else
        error( @sprintf("ERROR: unknown method = %s", method) )
    end

end

#=
Result from ABINIT: (30 Ry):
    Kinetic energy  =  2.16029548905902E-01
    Hartree energy  =  1.78287277495300E-01
    XC energy       = -2.76638068416253E-01
    Ewald energy    = -2.19685854008068E-02
    PspCore energy  = -3.96586964629223E-03
    Loc. psp. energy= -6.14021970150543E-01
    NL   psp  energy=  1.48929934899351E-01
    >>>>>>>>> Etotal= -3.73347732313342E-01

=#
