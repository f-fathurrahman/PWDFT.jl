mutable struct Hamiltonian{Txc<:AbstractXCCalculator}
    pw::PWGrid
    potentials::Potentials
    energies::Energies
    rhoe::Array{Float64,2} # spin dependent
    electrons::Electrons
    atoms::Atoms
    sym_info::SymmetryInfo
    rhoe_symmetrizer::RhoeSymmetrizer
    pspots::Array{PsPot_GTH,1}
    pspotNL::PsPotNL
    xcfunc::String
    xc_calc::Txc
    #xc_calc::Union{LibxcXCCalculator,XCCalculator}
    #xc_calc::AbstractXCCalculator
    ik::Int64   # current kpoint index
    ispin::Int64 # current spin index
end

"""
    Ham = Hamiltonian(...)

Create an instance of `Hamiltonian`.

Mandatory arguments:

- `atoms`: an instance of `Atoms`
- `pspfiles`: an array of `String` containing path to the pseudopotential files.
  They must be given such that they match the order of `atoms.SpeciesSymbols`.
- `ecutwfc`: energy cutoff for constructing plane wave basis set.

Thw following is the most commonly used optional arguments:
- `Nspin`: number of spin components. Only `Nspin=1` (no spin polarization) and
  `Nspin=2` (spin polarized) are supported.
- `meshk`: an array of three integers specifying grid size for Monkhorst-Pack grid.
  They will be reduced according the cell symmetry.
  This should be specified for crystalline systems.
- `shiftk`: specifying whether the grid will be shifted or not (using PWSCF convention).
  Presently only `[0,0,0]` is tested.
- `extra_states`: number of empty states. This must be specified for metallic systems.
- `xcfunc`: currently only support "VWN" or "PBE".
- `use_symmetry`: whether symmetry is used or not.
- `Ns_`: a tuple of three integers for overriding FFT grid.
- `use_xc_internal`: if false then Libxc will be used.
"""
function Hamiltonian( atoms::Atoms, pspfiles::Array{String,1},
                      ecutwfc::Float64 ;
                      Nspin=1,
                      meshk=[1,1,1], shiftk=[0,0,0], time_reversal=true,
                      Ns_=(0,0,0),
                      kpoints=nothing,
                      kpts_str="",
                      xcfunc="VWN",
                      use_xc_internal=false,
                      extra_states=0,
                      use_symmetry=true )

    if use_symmetry == false
        sym_info = SymmetryInfo()
    else
        sym_info = SymmetryInfo(atoms)
    end

    # kpoints
    if kpoints == nothing
        if kpts_str == ""
            # automatic generation of kpoints
            kpoints = KPoints( atoms, meshk, shiftk, sym_info.s, time_reversal=time_reversal )
        else
            # use the given kpoints
            kpoints = kpoints_from_string( atoms, kpts_str )
        end
    else
        @assert typeof(kpoints) == KPoints
    end

    # Initialize plane wave grids
    pw = PWGrid( ecutwfc, atoms.LatVecs, kpoints=kpoints, Ns_=Ns_ )

    Nspecies = atoms.Nspecies
    if Nspecies != size(pspfiles)[1]
        error( @sprintf("Length of pspfiles is not equal to %d\n", Nspecies) )
    end

    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume
    G2 = pw.gvec.G2
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r

    strf = calc_strfact( atoms, pw )

    #
    # Initialize pseudopotentials and local potentials
    #
    Pspots = Array{PsPot_GTH}(undef,Nspecies)
    G2_shells = pw.gvec.G2_shells
    Ngl = length(G2_shells)
    Vgl = zeros(Float64, Ngl)
    idx_g2shells = pw.gvec.idx_g2shells
    V_Ps_loc = zeros(Float64, Npoints)
    Vg = zeros(ComplexF64, Npoints)

    for isp = 1:Nspecies
        Pspots[isp] = PsPot_GTH( pspfiles[isp] )
        psp = Pspots[isp]
        #
        for igl = 1:Ngl
            Vgl[igl] = eval_Vloc_G( psp, G2_shells[igl] )
        end
        #
        for ig = 1:Ng
            ip = idx_g2r[ig]
            igl = idx_g2shells[ig]
            Vg[ip] = strf[ig,isp] * Vgl[igl]
        end
        #
        V_Ps_loc[:] = V_Ps_loc[:] + real( G_to_R(pw, Vg) ) * Npoints / CellVolume
    end

    # other potential terms are set to zero
    V_Hartree = zeros( Float64, Npoints )
    V_xc = zeros( Float64, Npoints, Nspin )
    V_loc_tot = zeros( Float64, Npoints, Nspin )
    potentials = Potentials( V_Ps_loc, V_Hartree, V_xc, V_loc_tot )
    #
    energies = Energies()
    #
    rhoe = zeros( Float64, Npoints, Nspin )

    electrons = Electrons( atoms, Pspots, Nspin=Nspin, Nkpt=kpoints.Nkpt,
                           Nstates_empty=extra_states )

    # NL pseudopotentials
    pspotNL = PsPotNL( atoms, pw, Pspots, check_norm=false )

    atoms.Zvals = get_Zvals( Pspots )

    ik = 1
    ispin = 1

    if sym_info.Nsyms > 1
        rhoe_symmetrizer = RhoeSymmetrizer( atoms, pw, sym_info )
    else
        rhoe_symmetrizer = RhoeSymmetrizer() # dummy rhoe_symmetrizer
    end

    if use_xc_internal
        xc_calc = XCCalculator()
    else
        # Using Libxc is the default
        if xcfunc == "SCAN"
            xc_calc = LibxcXCCalculator(is_metagga=true, Npoints=Npoints, Nspin=Nspin)
        else
            xc_calc = LibxcXCCalculator()
        end
    end

    return Hamiltonian( pw, potentials, energies, rhoe,
                        electrons, atoms, sym_info, rhoe_symmetrizer,
                        Pspots, pspotNL, xcfunc, xc_calc, ik, ispin )
end


#
# No pspfiles given. Use Coulomb potential (all electrons)
#
# WARNING: This is not tested extensively
function Hamiltonian( atoms::Atoms, ecutwfc::Float64;
                      Nspin=1,
                      meshk=[1,1,1],
                      shiftk=[0,0,0],
                      kpts_str="",
                      kpoints=nothing,   
                      xcfunc="VWN",
                      use_xc_internal=false,
                      extra_states=0 )

    sym_info = SymmetryInfo(atoms)

    # kpoints
    if kpoints == nothing
        if kpts_str == ""
            # automatic generation of kpoints
            kpoints = KPoints( atoms, meshk, shiftk )
        else
            # use the given kpoints
            kpoints = kpoints_from_string( atoms, kpts_str )
        end
    else
        @assert typeof(kpoints) == KPoints
    end


    # Initialize plane wave grids
    pw = PWGrid( ecutwfc, atoms.LatVecs, kpoints=kpoints )

    Nspecies = atoms.Nspecies

    Npoints = prod(pw.Ns)

    # dummy pspots
    Pspots = Array{PsPot_GTH}(undef,1)
    Pspots[1] = PsPot_GTH()

    #
    # Coulomb potential as local potential
    #
    Zvals = get_Zatoms(atoms)
    strf = calc_strfact( atoms, pw )
    V_Ps_loc = init_V_coulomb_G( pw, strf, Zvals )

    # other potentials are set to zero
    V_Hartree = zeros( Float64, Npoints )
    V_xc = zeros( Float64, Npoints, Nspin )
    V_loc_tot = zeros( Float64, Npoints, Nspin )
    potentials = Potentials( V_Ps_loc, V_Hartree, V_xc, V_loc_tot )
    #
    energies = Energies()

    Zatoms = get_Zatoms( atoms )
    atoms.Zvals = Zatoms
    
    # use Zatoms as Zvals
    electrons = Electrons( atoms, Zatoms, Nspin=Nspin, Nkpt=kpoints.Nkpt,
                           Nstates_empty=extra_states )

    rhoe = zeros(Float64,Npoints,Nspin)

    pspotNL = PsPotNL()
    
    ik = 1
    ispin = 1

    if sym_info.Nsyms > 1
        rhoe_symmetrizer = RhoeSymmetrizer( atoms, pw, sym_info )
    else
        rhoe_symmetrizer = RhoeSymmetrizer() # dummy rhoe_symmetrizer
    end

    if use_xc_internal
        xc_calc = XCCalculator()
    else
        # Using Libxc is the default
        xc_calc = LibxcXCCalculator()
    end

    return Hamiltonian( pw, potentials, energies, rhoe,
                        electrons, atoms, sym_info, rhoe_symmetrizer,
                        Pspots, pspotNL, xcfunc, xc_calc, ik, ispin )
end

# FIXME: should be merged with Array{Float64,2} version.
# psiks is needed for metagga
function update!( Ham::Hamiltonian, psiks::BlochWavefunc, rhoe::Array{Float64,1} )
    # assumption Nspin = 1
    Ham.rhoe[:,1] = rhoe
    Vxc = zeros(size(rhoe,1))
    Ham.potentials.Hartree = real( G_to_R( Ham.pw, Poisson_solve(Ham.pw, rhoe) ) )    
    if Ham.xcfunc == "SCAN"
        calc_Vxc_SCAN!( Ham.xc_calc, Ham.pw, psiks, rhoe, Vxc )
        Ham.potentials.XC[:,1] = Vxc[:]
    elseif Ham.xcfunc == "PBE"
        Ham.potentials.XC[:,1] = calc_Vxc_PBE( Ham.xc_calc, Ham.pw, rhoe )
    else  # VWN is the default
        Ham.potentials.XC[:,1] = calc_Vxc_VWN( Ham.xc_calc, rhoe )
    end
    
    Npoints = prod(Ham.pw.Ns)
    for ip = 1:Npoints
        Ham.potentials.Total[ip,1] = Ham.potentials.Ps_loc[ip] + Ham.potentials.Hartree[ip] +
                                     Ham.potentials.XC[ip,1]  
    end
    return
end


"""
    update!(Ham, psiks, rhoe)

Update Ham.rhoe and calculate Hartree and XC potentials for given `rhoe` in real space.
"""
function update!(Ham::Hamiltonian, psiks::BlochWavefunc, rhoe::Array{Float64,2})
    Nspin = size(rhoe)[2]
    if Nspin == 1
        update!(Ham, psiks, rhoe[:,1])
        return
    end
    # FIXME: spinpol MetaGGA is not yet implemented
    Ham.rhoe = rhoe[:,:]
    Rhoe_total = Ham.rhoe[:,1] + Ham.rhoe[:,2] # Nspin is 2
    Ham.potentials.Hartree = real( G_to_R( Ham.pw, Poisson_solve(Ham.pw, Rhoe_total) ) )
    if Ham.xcfunc == "PBE"
        Ham.potentials.XC = calc_Vxc_PBE( Ham.xc_calc, Ham.pw, rhoe )
    else  # VWN is the default
        Ham.potentials.XC = calc_Vxc_VWN( Ham.xc_calc, rhoe )
    end
    Npoints = prod(Ham.pw.Ns)
    for ispin = 1:Nspin
        for ip = 1:Npoints
            Ham.potentials.Total[ip,ispin] = Ham.potentials.Ps_loc[ip] + Ham.potentials.Hartree[ip] +
                                             Ham.potentials.XC[ip,ispin]  
        end
    end
    return
end

# For compatibility
function update!( Ham::Hamiltonian, rhoe::Array{Float64,1} )
    # assumption Nspin = 1
    Ham.rhoe[:,1] = rhoe
    Vxc = zeros(size(rhoe,1))
    Ham.potentials.Hartree = real( G_to_R( Ham.pw, Poisson_solve(Ham.pw, rhoe) ) )    
    if Ham.xcfunc == "PBE"
        Ham.potentials.XC[:,1] = calc_Vxc_PBE( Ham.xc_calc, Ham.pw, rhoe )
    else  # VWN is the default
        Ham.potentials.XC[:,1] = calc_Vxc_VWN( Ham.xc_calc, rhoe )
    end
    Npoints = prod(Ham.pw.Ns)
    for ip = 1:Npoints
        Ham.potentials.Total[ip,1] = Ham.potentials.Ps_loc[ip] + Ham.potentials.Hartree[ip] +
                                     Ham.potentials.XC[ip,1]  
    end
    return
end

function update!(Ham::Hamiltonian, rhoe::Array{Float64,2})
    Nspin = size(rhoe)[2]
    if Nspin == 1
        update!(Ham, rhoe[:,1])
        return
    end
    Ham.rhoe = rhoe[:,:]
    Rhoe_total = Ham.rhoe[:,1] + Ham.rhoe[:,2] # Nspin is 2
    Ham.potentials.Hartree = real( G_to_R( Ham.pw, Poisson_solve(Ham.pw, Rhoe_total) ) )
    if Ham.xcfunc == "PBE"
        Ham.potentials.XC = calc_Vxc_PBE( Ham.xc_calc, Ham.pw, rhoe )
    else  # VWN is the default
        Ham.potentials.XC = calc_Vxc_VWN( Ham.xc_calc, rhoe )
    end
    Npoints = prod(Ham.pw.Ns)
    for ispin = 1:Nspin
        for ip = 1:Npoints
            Ham.potentials.Total[ip,ispin] = Ham.potentials.Ps_loc[ip] + Ham.potentials.Hartree[ip] +
                                             Ham.potentials.XC[ip,ispin]  
        end
    end
    return
end


include("Hamiltonian_io.jl")