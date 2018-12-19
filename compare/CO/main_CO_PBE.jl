function main( ; method="SCF" )

    # Atoms
    atoms = Atoms( xyz_file="../structures/CO.xyz",
                   LatVecs = gen_lattice_cubic(16.0) )

    # Initialize Hamiltonian
    ecutwfc = 15.0
    pspfiles = ["../pseudopotentials/pbe_gth/C-q4.gth",
                "../pseudopotentials/pbe_gth/O-q6.gth"]
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, xcfunc="PBE" )
    println(Ham)

    if method == "SCF"
        KS_solve_SCF!( Ham, mix_method="anderson" )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham, verbose=true )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    else
        error( @sprintf("ERROR: unknown method = %s", method) )
    end

end

#=
    Kinetic energy  =  1.27516710648293E+01
    Hartree energy  =  1.88243902547401E+01
    XC energy       = -4.85795922450072E+00
    Ewald energy    =  2.26310215311565E+00
    PspCore energy  =  2.21150310652737E-04
    Loc. psp. energy= -5.21312785756274E+01
    NL   psp  energy=  2.56490073237464E+00
    >>>>>>>>> Etotal= -2.05849524447578E+01


Kinetic    energy:      12.6571089076
Ps_loc     energy:     -51.9688274308
Ps_nloc    energy:       2.5943370534
Hartree    energy:      18.7034268280
XC         energy:      -4.8329259033
-------------------------------------
Electronic energy:     -22.8468805451
NN         energy:       2.2631029868
-------------------------------------
Total      energy:     -20.5837775583


Kinetic    energy:      12.6573215345
Ps_loc     energy:     -51.9690474160
Ps_nloc    energy:       2.5947155672
Hartree    energy:      18.7028490834
XC         energy:      -4.8328186269
-------------------------------------
Electronic energy:     -22.8469798577
NN         energy:       2.2631029868
-------------------------------------
Total      energy:     -20.5838768709
=#
