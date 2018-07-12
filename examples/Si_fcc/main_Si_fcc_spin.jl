using Printf
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
    pspfiles = ["../../pseudopotentials/pade_gth/Si-q4.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, xcfunc="VWN",
                         meshk=[3,3,3], verbose=true, Nspin=2, extra_states=0 )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

    #
    # Solve the KS problem
    #
    #KS_solve_Emin_PCG!( Ham )
    KS_solve_SCF!( Ham, mix_method="anderson" )
    #KS_solve_SCF_smearing!( Ham, update_psi="LOBPCG", mix_method="simple" )

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

#=
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

--------------

With smearing 0.01 Ha
Total energy components

Kinetic    energy:       3.2023421562
Ps_loc     energy:      -2.2762398006
Ps_nloc    energy:       1.5915987675
Hartree    energy:       0.5322060942
XC         energy:      -2.2975477892
-------------------------------------
Electronic energy:       0.7523594282
NN         energy:      -8.3979274007
-------------------------------------
Total      energy:      -7.6455679725

------------

ABINIT

Kinetic energy  =  3.18425700948638E+00
Hartree energy  =  5.65422931914669E-01
XC energy       = -2.42294101779165E+00
Ewald energy    = -8.39792740071415E+00
PspCore energy  = -2.09890966765802E-01
Loc. psp. energy= -2.12992717173527E+00
NL   psp  energy=  1.55983881823907E+00
>>>>> Internal E= -7.85116779736675E+00

-kT*entropy     = -4.52029451264488E-03
>>>>>>>>> Etotal= -7.85568809187940E+00

PWSCF

!    total energy              =     -15.70835819 Ry = -7.854179095
     Harris-Foulkes estimate   =     -15.70835757 Ry
     estimated scf accuracy    <       0.00000013 Ry
     one-electron contribution =       4.79575036 Ry
     hartree contribution      =       1.14632766 Ry
     xc contribution           =      -4.85295131 Ry
     ewald contribution        =     -16.79585054 Ry
     smearing contrib. (-TS)   =      -0.00163437 Ry

=#
