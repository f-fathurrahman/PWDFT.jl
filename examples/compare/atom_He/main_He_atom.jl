function main( ; method="SCF" )

    # Atoms
    atoms = Atoms(xyz_file=joinpath(DIR_STRUCTURES, "He.xyz"),
                  LatVecs=gen_lattice_sc(16.0))

    # Initialize Hamiltonian
    pspfiles = [joinpath(DIR_PSP, "He-q2.gth")]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc )

    if method == "SCF"
        KS_solve_SCF!( Ham )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    else
        error( @sprintf("Unknown method = %s", method) )
    end

end

#=
ABINIT (result):
    Kinetic energy  =  2.18802875961751E+00
    Hartree energy  =  1.50108932871126E+00
    XC energy       = -9.02904237573505E-01
    Ewald energy    = -3.54662184935077E-01
    PspCore energy  = -1.69161860064652E-06
    Loc. psp. energy= -5.08088416764117E+00
    NL   psp  energy=  0.00000000000000E+00
    >>>>>>>>> Etotal= -2.64933419343958E+00

Total energy components
Kinetic    energy:       2.1883489615
Ps_loc     energy:      -5.0813492117
Ps_nloc    energy:       0.0000000000
Hartree    energy:       1.5013376586
XC         energy:      -0.9038956763
-------------------------------------
Electronic energy:      -2.2955582679
NN         energy:      -0.3546621849
-------------------------------------
Total      energy:      -2.6502204529

=#
