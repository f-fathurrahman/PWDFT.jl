# A helper function to print forces
function _print_forces(atoms::Atoms, F, title_str)
    Natoms = atoms.Natoms
    atsymbs = atoms.atsymbs
    println()
    println(title_str)
    for ia in 1:Natoms
        @printf("%2s %18.10f %18.10f %18.10f\n", atsymbs[ia], F[1,ia], F[2,ia], F[3,ia])
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

    if isnothing(filename)
        export_to_script(pwinput)
    else
        export_to_script(pwinput, filename="script_"*filename*".jl")
    end

    if pwinput.occupations == "smearing"
        Ham.electrons.use_smearing = true
        Ham.electrons.kT = pwinput.degauss*0.5 # convert from Ry to Ha
    end

    if pwinput.nspin == 2
        starting_magn = pwinput.starting_magnetization
    else
        starting_magn = nothing
    end

    #electrons_scf!(Ham, psiks, NiterMax=100, use_smearing=use_smearing, kT=kT, betamix=0.1)
    
    psiks = rand_BlochWavefunc(Ham)
    electrons_scf_G!(
        Ham,
        psiks=psiks,
        NiterMax=100,
        betamix=0.1,
        starting_magn=starting_magn
    )
    # XXX "Restarting" from previous calculation can be done by specifying Rhoe=Ham.rhoe

    
    #KS_solve_SCF!(Ham, psiks, use_smearing=use_smearing, kT=kT, betamix=0.1)
    # Not yet working for smearing
    #KS_solve_Emin_PCG!(Ham, psiks, NiterMax=100)

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


# This is a left-over from testing electrons_Emin_Haux!
# Some care is given for GTH/HGH pseudopotentials in case of magnetic calculation.
function prepare_Ham_from_pwinput(; filename=nothing)

    Ham, pwinput = init_Ham_from_pwinput(filename=filename);

    # Compute this once and for all
    Ham.energies.NN = calc_E_NN(Ham.atoms);
    
    # We need to set some parameters manually:
    use_smearing = false
    kT = 0.0
    if pwinput.occupations == "smearing"
        use_smearing = true
        kT = pwinput.degauss*0.5 # convert from Ry to Ha
        Ham.electrons.use_smearing = true
        Ham.electrons.kT = kT
    end

    if pwinput.nspin == 2
        starting_magn = pwinput.starting_magnetization
    else
        starting_magn = nothing
    end

    Nspin = Ham.electrons.Nspin_channel
    # Initialize electronic variables: `psiks` and `Haux`:
    Random.seed!(1234)
    psiks = rand_BlochWavefunc(Ham)
    # XXX: this is not really needed, should be able to pass nothing to

    # Initial Rhoe, taking into account starting magnetization
    #
    # XXX: This is needed for GTH/HGH pspots because rho_atom is not available
    # some heuristics
    is_using_gth_analytic = (eltype(Ham.pspots) == PsPot_GTH)
    is_using_gth_numeric = contains(uppercase(Ham.pspots[1].pspfile), "GTH")
    is_using_hgh_numeric = contains(uppercase(Ham.pspots[1].pspfile), "HGH")
    no_atomic_rhoe = is_using_gth_analytic || is_using_gth_numeric || is_using_hgh_numeric
    if Nspin == 2 && no_atomic_rhoe
        _, _ = _prepare_scf!(Ham, psiks)
        _, _ = update_from_rhoe!( Ham, psiks, Ham.rhoe )
        Rhoe_tot = Ham.rhoe[:,1] + Ham.rhoe[:,2]
        magn = starting_magn .* ones(size(Rhoe_tot)) / Ham.pw.CellVolume
        Ham.rhoe[:,1] .= 0.5*(Rhoe_tot + magn)
        Ham.rhoe[:,2] .= 0.5*(Rhoe_tot - magn)
        # Update again the Hamiltonian
        _, _ = update_from_rhoe!( Ham, psiks, Ham.rhoe )
    else
        # for everything else atomic rhoe from pspots should be available
        _, _ = _prepare_scf!(Ham, psiks, starting_magn=starting_magn)
    end

    return Ham
end
