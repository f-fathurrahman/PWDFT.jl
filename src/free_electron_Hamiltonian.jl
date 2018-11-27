function free_electron_Hamiltonian(
        atoms::Atoms, pspfiles::Array{String,1},
        ecutwfc::Float64 ;
        Nspin = 1,
        meshk = [1,1,1], shiftk = [0,0,0],
        kpoints = nothing,
        xcfunc = "VWN",
        verbose = false,
        extra_states = 0
    )

    # kpoints
    if kpoints == nothing
        kpoints = KPoints( atoms, meshk, shiftk )
    else
        @assert typeof(kpoints) == KPoints
    end

    # Initialize plane wave grids
    pw = PWGrid( ecutwfc, atoms.LatVecs, kpoints=kpoints )
    if verbose
        println(pw)
        println(kpoints)
    end

    Nspecies = atoms.Nspecies
    if Nspecies != size(pspfiles)[1]
        error( @sprintf("Length of pspfiles is not equal to %d\n", Nspecies) )
    end

    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume
    G2 = pw.gvec.G2
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r

    #
    # Initialize pseudopotentials and local potentials
    #
    Vg = zeros(ComplexF64, Npoints)
    V_Ps_loc = zeros(Float64, Npoints)

    Pspots = Array{PsPot_GTH}(undef,Nspecies)

    for isp = 1:Nspecies
        Pspots[isp] = PsPot_GTH( pspfiles[isp] )
    end

    # other potential terms are set to zero
    V_Hartree = zeros( Float64, Npoints )
    V_XC = zeros( Float64, Npoints, Nspin )
    potentials = Potentials( V_Ps_loc, V_Hartree, V_XC )
    #
    energies = Energies()
    #
    rhoe = zeros( Float64, Npoints, Nspin )

    electrons = Electrons( atoms, Pspots, Nspin=Nspin, Nkpt=kpoints.Nkpt,
                           Nstates_empty=extra_states )
    if verbose
        println(electrons)
    end

    # NL pseudopotentials
    pspotNL = PsPotNL( pw, atoms, Pspots, kpoints, check_norm=false )

    atoms.Zvals = get_Zvals( Pspots )

    ik = 1
    ispin = 1
    return Hamiltonian( pw, potentials, energies, rhoe,
                        electrons, atoms, Pspots, pspotNL, xcfunc, ik, ispin )
end