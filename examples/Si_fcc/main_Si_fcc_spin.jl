using PWDFT

function test_main()
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(10.2631)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)
    write_xsf( "TEMP_Si.xsf", atoms )

    # Initialize Hamiltonian
    #pspfiles = ["../../pseudopotentials/pade_gth/Si-q4.gth"]
    pspfiles = ["../../pseudopotentials/pbe_gth/Si-q4.gth"]
    ecutwfc_Ry = 30.0
    Ham = PWHamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, xcfunc="PBE",
                         meshk=[3,3,3], verbose=true, Nspin=2 )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

    #
    # Solve the KS problem
    #
    #KS_solve_SCF!( Ham, update_psi="LOBPCG" )
    KS_solve_SCF_smearing!( Ham, update_psi="LOBPCG" )

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

@time test_main()

"""

No smearing:

Total energy components

Kinetic    energy:       3.1744863175
Ps_loc     energy:      -2.3423254958
Ps_nloc    energy:       1.5694366997
Hartree    energy:       0.5730118848
XC         energy:      -2.4296659652
-------------------------------------
Electronic energy:       0.5449434409
NN         energy:      -8.3979274007
-------------------------------------
Total      energy:      -7.8529839598
"""
