function main( ; method="SCF" )

    # Atoms
    atoms = Atoms( xyz_file="../structures/H2.xyz",
                   LatVecs=gen_lattice_sc(16.0) )

    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth"]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc,
                         Nspin=2, extra_states=0 )
    println(Ham)

    if method == "SCF"
        KS_solve_SCF!( Ham, mix_method="anderson" )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    else
        error( @sprintf("ERROR: unknown method = %s", method) )
    end

end

#=
For 30 Ry

ABINIT result:
Kinetic energy  =  1.01004063629599E+00
Hartree energy  =  9.01545108206733E-01
XC energy       = -6.31436411843721E-01
Ewald energy    =  3.13170052325859E-01
PspCore energy  = -1.26742500464741E-06
Loc. psp. energy= -2.71283251550643E+00
NL   psp  energy=  0.00000000000000E+00
>>>>> Internal E= -1.11951439794657E+00
-kT*entropy     = -7.09088279385096E-19
>>>>>>>>> Etotal= -1.11951439794657E+00

Kinetic    energy:       1.0100083461
Ps_loc     energy:      -2.7127890544
Ps_nloc    energy:       0.0000000000
Hartree    energy:       0.9015208969
XC         energy:      -0.6314252227
-------------------------------------
Electronic energy:      -1.4326850340
NN         energy:       0.3131700523
-------------------------------------
Total      energy:      -1.1195149817
=#

