function main( ; method="SCF" )

    # Atoms
    atoms = Atoms( xyz_file="../structures/CO.xyz",
                   LatVecs = gen_lattice_sc(16.0) )

    # Initialize Hamiltonian
    ecutwfc_Ry = 30.0
    pspfiles = ["../pseudopotentials/pade_gth/C-q4.gth",
                "../pseudopotentials/pade_gth/O-q6.gth"]
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5 )
    println(Ham)

    if method == "SCF"
        KS_solve_SCF!( Ham, mix_method="anderson" )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    elseif method == "TRDCM"
        KS_solve_TRDCM!( Ham, NiterMax=15 )

    else
        error( @sprintf("ERROR: unknown method = %s", method) )
    end

end

#=
    Kinetic energy  =  1.27126195273501E+01
    Hartree energy  =  1.87563651375840E+01
    XC energy       = -4.74110355787942E+00
    Ewald energy    =  2.26310215311565E+00
    PspCore energy  = -2.54693990241102E-04
    Loc. psp. energy= -5.21226082214705E+01
    NL   psp  energy=  2.59183822442630E+00
    >>>>>>>>> Etotal= -2.05400414308641E+01

Kinetic    energy:      12.7122820382
Ps_loc     energy:     -52.1220967254
Ps_nloc    energy:       2.5918770838
Hartree    energy:      18.7558153675
XC         energy:      -4.7410160310
-------------------------------------
Electronic energy:     -22.8031382669
NN         energy:       2.2631029868
-------------------------------------
Total      energy:     -20.5400352801
=#
