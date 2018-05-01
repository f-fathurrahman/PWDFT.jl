using PWDFT

function test_main( ; method="SCF" )
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
    write_xsf( "TEMP_Si.xsf", atoms )

    # Initialize Hamiltonian
    pspfiles = ["../../pseudopotentials/pade_gth/Si-q4.gth"]
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

#@time test_main(method="Emin")
@time test_main(method="SCF")
#@time test_main(method="DCM")


"""
FFT grid = (27,27,27)
!    total energy              =     -15.82591729 Ry = -7.912958645
     one-electron contribution =       4.64239409 Ry
     hartree contribution      =       1.15137054 Ry
     xc contribution           =      -4.82383139 Ry
     ewald contribution        =     -16.79585054 Ry = -8.39792527

FFT grid = (27,27,27)
ABINIT result
Kinetic energy  =  3.24281107903430E+00
Hartree energy  =  5.65725361752616E-01
XC energy       = -2.42086670957799E+00
Ewald energy    = -8.46648022654903E+00
PspCore energy  = -3.01899831461368E-01
Loc. psp. energy= -2.11915547489918E+00
NL   psp  energy=  1.59139219889343E+00
>>>>>>>>> Etotal= -7.90847360280721E+00

FFT grid = (27,27,27)
Kinetic    energy:       3.2107925351
Ps_loc     energy:      -2.4706316946
Ps_nloc    energy:       1.5806564220
Hartree    energy:       0.5830491953
XC         energy:      -2.4157133347
-------------------------------------
Electronic energy:       0.4881531232
NN         energy:      -8.3979258900
-------------------------------------
Total      energy:      -7.9097727668
"""