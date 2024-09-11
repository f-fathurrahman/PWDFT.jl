# A helper function to print forces
function _print_forces(atoms::Atoms, F, title_str)
    Natoms = atoms.Natoms
    atsymbs = atoms.atsymbs
    println()
    println(title_str)
    for ia in 1:Natoms
        @printf("%s %18.10f %18.10f %18.10f\n", atsymbs[ia], F[1,ia], F[2,ia], F[3,ia])
    end
    return
end

"""
This function roughly emulates pw.x binary in Quantum Espresso.
It reads pw.x input file given in keyword argument `filename`.
If no argument is given the default filename is `PWINPUT` which should be
present in the current directory.
"""
function my_pwx(; filename=nothing, do_export_data=false)
    Ham, pwinput = init_Ham_from_pwinput(filename=filename)

    write_xsf("ATOMS_from_pwinput.xsf", Ham.atoms)
    println(Ham)

    export_to_script(pwinput, filename="script_"*filename*".jl")

    # This will take into account whether the overlap operator is needed or not
    psiks = rand_BlochWavefunc(Ham)

    use_smearing = false
    kT = 0.0
    if pwinput.occupations == "smearing"
        use_smearing = true
        kT = pwinput.degauss*0.5 # convert from Ry to Ha
    end

    electrons_scf!(Ham, psiks, NiterMax=100, use_smearing=use_smearing, kT=kT, betamix=0.1)

    # Calculate forces: we want to display them by different contributions
    # in Ry/bohr to facilitate easier comparison with QE
    #
    # Factor of 2 to convert to Ry/bohr
    F_NN = 2*calc_forces_NN( Ham.pw, Ham.atoms )
    F_Ps_loc = 2*calc_forces_Ps_loc( Ham )
    F_Ps_nloc = 2*calc_forces_Ps_nloc( Ham, psiks )
    symmetrize_vector!( Ham.pw.LatVecs, Ham.sym_info, F_Ps_nloc )
    #
    F_nlcc = 2*calc_forces_nlcc(Ham)
    F_scf_corr = 2*calc_forces_scf_corr(Ham)
    #
    F_tot = F_NN + F_Ps_loc + F_Ps_nloc + F_nlcc + F_scf_corr

    _print_forces(Ham.atoms, F_tot, "Atomic forces (in Ry/bohr)")
    _print_forces(Ham.atoms, F_Ps_nloc, "Nonlocal pspot contribution to forces (in Ry/bohr)")
    _print_forces(Ham.atoms, F_NN, "Ionic contribution to forces (in Ry/bohr)")
    _print_forces(Ham.atoms, F_Ps_loc, "Local contribution to forces (in Ry/bohr)")
    _print_forces(Ham.atoms, F_nlcc, "Core correction contribution to forces (in Ry/bohr)")
    _print_forces(Ham.atoms, F_scf_corr, "SCF correction to forces (in Ry/bohr)")

    if do_export_data
        Serialization.serialize("Hamiltonian.dat", Ham)
        Serialization.serialize("psiks.dat", psiks)
        @info("")
        @info("Hamiltonian and psiks are serialized to files:")
        @info("Hamiltonian.dat and psiks.dat")
        @info("")
        @info("!!! Beware that FFTW plans are C-pointers.")
        @info("!!! You should not use them from serialized Hamiltonian")
    end

#=
    # Not yet working for smearing
    KS_solve_Emin_PCG!(Ham, psiks, NiterMax=100)
    energies = Ham.energies
    println()
    println(">>>> Final result:")
    println()
    
    println("-------------------------------------")
    println("Energy components in Ry")
    println("-------------------------------------")
    
    @printf("Kinetic    energy: %18.10f Ry\n", 2*energies.Kinetic )
    @printf("Ps_loc     energy: %18.10f Ry\n", 2*energies.Ps_loc )
    @printf("Ps_nloc    energy: %18.10f Ry\n", 2*energies.Ps_nloc )
    @printf("Hartree    energy: %18.10f Ry\n", 2*energies.Hartree )
    @printf("XC         energy: %18.10f Ry\n", 2*energies.XC )
    @printf("-TS              : %18.10f Ry\n", 2*energies.mTS)
    @printf("-------------------------------------\n")
    
    E_elec = energies.Kinetic + energies.Ps_loc + energies.Ps_nloc +
             energies.Hartree + energies.XC + energies.mTS
    
    @printf("Electronic energy: %18.10f Ry\n", 2*E_elec)
    @printf("NN         energy: %18.10f Ry\n", 2*energies.NN )
    @printf("-------------------------------------\n")

    E_total = E_elec + energies.NN
    @printf("! Total = %18.10f Ry\n", 2*E_total)
=#
    return

end
