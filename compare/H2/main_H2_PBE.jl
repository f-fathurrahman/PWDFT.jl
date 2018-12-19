function main( ; method="SCF" )

    # Atoms
    atoms = Atoms( xyz_file="../structures/H2.xyz",
                   LatVecs=gen_lattice_sc(16.0) )

    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pbe_gth/H-q1.gth"]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, xcfunc="PBE" )
    println(Ham)

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
    Kinetic energy  =  1.03846627435554E+00
    Hartree energy  =  9.16452467490984E-01
    XC energy       = -6.66546679987766E-01
    Ewald energy    =  3.13170052325859E-01
    PspCore energy  = -1.32876358319460E-06
    Loc. psp. energy= -2.74851348731632E+00
    NL   psp  energy=  0.00000000000000E+00
    >>>>>>>>> Etotal= -1.14697270189529E+00

Kinetic    energy:       1.0022039947
Ps_loc     energy:      -2.7005240671
Ps_nloc    energy:       0.0000000000
Hartree    energy:       0.8935857782
XC         energy:      -0.6547008647
-------------------------------------
Electronic energy:      -1.4594351588
NN         energy:       0.3131700523
-------------------------------------
Total      energy:      -1.1462651065
=#

