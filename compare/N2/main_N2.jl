function main( ; method="SCF" )

    # Atoms
    atoms = Atoms( xyz_file="../structures/N2.xyz",
                   LatVecs=gen_lattice_sc(16.0) )

    # Initialize Hamiltonian
    ecutwfc_Ry = 30.0
    pspfiles = ["../pseudopotentials/pade_gth/N-q5.gth"]
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5 )
    println(Ham)

    if method == "SCF"
        KS_solve_SCF!( Ham, mix_method="rpulay", betamix=0.2 )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    else
        error( @sprintf("ERROR: unknown method = %s", method) )
    end

end

#=
Kinetic    energy:      11.5329676573
Ps_loc     energy:     -47.3445263895
Ps_nloc    energy:       2.3268197661
Hartree    energy:      17.1133647701
XC         energy:      -4.5072691732
-------------------------------------
Electronic energy:     -20.8786433691
NN         energy:       1.7906353998
-------------------------------------
Total      energy:     -19.0880079692

    Kinetic energy  =  1.15327757124958E+01
    Hartree energy  =  1.71131691439872E+01
    XC energy       = -4.50723738034607E+00
    Ewald energy    =  1.79063539962452E+00
    PspCore energy  = -7.02139897582120E-05
    Loc. psp. energy= -4.73441139820555E+01
    NL   psp  energy=  2.32683843351565E+00
    >>>>>>>>> Etotal= -1.90880028867682E+01


PWSCF result for 30 Ry
!    total energy              =     -38.17145300 Ry = -19.0857265 Ha
     one-electron contribution =     -66.96020106 Ry
     hartree contribution      =      34.21524310 Ry
     xc contribution           =      -9.00776423 Ry
     ewald contribution        =       3.58126919 Ry =   1.790634595 Ha
=#
