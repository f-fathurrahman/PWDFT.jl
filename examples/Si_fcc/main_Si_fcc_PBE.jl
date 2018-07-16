function main( ; method="SCF" )
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(5.431*ANG2BOHR)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)

    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pbe_gth/Si-q4.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5,
                         xcfunc="PBE", meshk=[3,3,3], verbose=true )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

    #
    # Solve the KS problem
    #
    if method == "SCF"
        KS_solve_SCF!( Ham, update_psi="PCG" )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham, verbose=true )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    else
        println("ERROR: unknown method = ", method)
    end

    println("\nTotal energy components")
    println(Ham.energies)

end

#=
Kinetic energy  =  3.19227108452527E+00
Hartree energy  =  5.74882034691866E-01
XC energy       = -2.42725552640522E+00
Ewald energy    = -8.39792740071415E+00
PspCore energy  = -2.09890966765802E-01
Loc. psp. energy= -2.14639614323146E+00
NL   psp  energy=  1.56065513075350E+00
>>>>> Internal E= -7.85366178714600E+00
-kT*entropy     = -3.43186293502298E-04
>>>>>>>>> Etotal= -7.85400497343950E+00

Kinetic energy  =  3.19382283897098E+00
Hartree energy  =  5.76523077574832E-01
XC energy       = -2.42800413076515E+00
Ewald energy    = -8.39792740071415E+00
PspCore energy  = -2.09890966765802E-01
Loc. psp. energy= -2.14931165868622E+00
NL   psp  energy=  1.56084552734338E+00
>>>>>>>>> Etotal= -7.85394271304213E+00

Total energy components (nospin)
Kinetic    energy:       3.1725477753
Ps_loc     energy:      -2.3381319811
Ps_nloc    energy:       1.5676041118
Hartree    energy:       0.5715793505
XC         energy:      -2.4291537146
-------------------------------------
Electronic energy:       0.5444455419
NN         energy:      -8.3979258900
-------------------------------------
Total      energy:      -7.8534803481
=#
