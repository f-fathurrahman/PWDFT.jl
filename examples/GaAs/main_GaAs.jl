using PWDFT

function test_main( ; method="SCF" )
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        2

        Ga  0.0   0.0   0.0
        As  0.25  0.25  0.25
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(5.6537*ANG2BOHR)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)
    write_xsf( "TEMP_GaAs.xsf", atoms )

    # Initialize Hamiltonian
    pspfiles = ["../../pseudopotentials/pade_gth/Ga-q3.gth",
                "../../pseudopotentials/pade_gth/As-q5.gth"]
    ecutwfc_Ry = 30.0
    Ham = PWHamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, meshk=[3,3,3], verbose=true )

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

    Nstates = Ham.electrons.Nstates
    ebands = Ham.electrons.ebands
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    k = Ham.pw.gvecw.kpoints.k
    
    println("\nBand energies:")
    for ik = 1:Nkpt
        @printf("%d k = [%f,%f,%f]\n", ik, k[1,ik], k[2,ik], k[3,ik])
        for ist = 1:Nstates
            @printf("%8d  %18.10f = %18.10f eV\n", ist, ebands[ist,ik], ebands[ist,ik]*Ry2eV*2)
        end
    end
    
    println("\nTotal energy components")
    println(Ham.energies)

end

@time test_main(method="Emin")
@time test_main(method="SCF")
#@time test_main(method="DCM")

"""
Kinetic energy  =  3.26604416655604E+00
Hartree energy  =  8.16569937191940E-01
XC energy       = -2.39970840396867E+00
Ewald energy    = -8.42120578386954E+00
PspCore energy  =  3.78008143191317E-01
Loc. psp. energy= -3.13594871632497E+00
NL   psp  energy=  8.58003632435184E-01
>>>>>>>>> Etotal= -8.63823702478870E+00

PCG result
Kinetic    energy:       3.2667303969
Ps_loc     energy:      -2.7591418239
Ps_nloc    energy:       0.8579762330
Hartree    energy:       0.8240957924
XC         energy:      -2.4078020740
-------------------------------------
Electronic energy:      -0.2181414755
NN         energy:      -8.4212062785
-------------------------------------
Total      energy:      -8.6393477541

SCF result
Kinetic    energy:       3.2667262728
Ps_loc     energy:      -2.7592009319
Ps_nloc    energy:       0.8580163543
Hartree    energy:       0.8241183786
XC         energy:      -2.4078096283
-------------------------------------
Electronic energy:      -0.2181495546
NN         energy:      -8.4212062785
-------------------------------------
Total      energy:      -8.6393558331
"""
