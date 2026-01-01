include("HamiltonianOptions.jl")

mutable struct Hamiltonian{Tpsp<:AbstractPsPot}
    pw::PWGrid
    potentials::Potentials
    energies::Energies
    rhoe::Array{Float64,2} # spin dependent
    rhoe_core::Union{Nothing,Array{Float64,2}}
    electrons::Electrons
    atoms::Atoms
    sym_info::SymmetryInfo
    rhoe_symmetrizer::RhoeSymmetrizer
    pspots::Vector{Tpsp}
    pspotNL::Union{PsPotNL,PsPotNL_UPF} # TODO: rename or add alias: PsPotNL to PsPotNL_GTH
    xcfunc::String
    xc_calc::Union{LibxcXCCalculator,XCCalculator}
    ik::Int64   # current kpoint index
    ispin::Int64 # current spin index
    need_overlap::Bool
    options::HamiltonianOptions # copy
end


function Hamiltonian(
    atoms::Atoms,
    pspots::Vector{Tpsp},
    ecutwfc::Float64,
    options::HamiltonianOptions
) where Tpsp <: AbstractPsPot

    domag = false # this is only have meaning in case of lspinorb (with noncollinear)
    if options.noncollinear
        if !isnothing(options.starting_magn)
            domag = true
        end
    end

    if options.use_symmetry == false
        sym_info = SymmetryInfo()
    else
        magnetic_sym = options.noncollinear && domag
        #time_reversal = .NOT. noinv .AND. .NOT. magnetic_sym
        if options.noncollinear
            s_mag = options.starting_magn
            angle1 = options.angle1
            angle2 = options.angle2
            m_loc = zeros(Float64, 3, atoms.Natoms)
            for ia in 1:atoms.Natoms
                isp = atoms.atm2species[ia]
                m_loc[1,ia] = s_mag[isp] * sind(angle1[isp]) * cosd(angle2[isp])
                m_loc[2,ia] = s_mag[isp] * sind(angle1[isp]) * sind(angle2[isp])
                m_loc[3,ia] = s_mag[isp] * cosd(angle1[isp])
            end
        else
            m_loc = nothing
        end
        sym_info = SymmetryInfo(atoms, magnetic_sym = magnetic_sym, m_loc = m_loc)
    end

    @assert options.dual >= 4.0

    # kpoints
    if isnothing(options.kpoints)
        if options.kpts_str == ""
            # automatic generation of kpoints
            kpoints = KPoints( atoms,
                options.meshk, options.shiftk,
                sym_info.s, time_reversal=options.time_reversal )
        else
            # use the given kpoints
            kpoints = kpoints_from_string( atoms, options.kpts_str )
        end
    else
        @assert typeof(kpoints) == KPoints
    end

    # Initialize plane wave grids
    pw = PWGrid( ecutwfc,
        atoms.LatVecs, kpoints=kpoints,
        Ns_=options.Ns,
        dual=options.dual
    )

    Nspecies = atoms.Nspecies
    @assert size(pspots,1) == Nspecies

    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume
    G2 = pw.gvec.G2
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r

    strf = calc_strfact( atoms, pw )

    #
    # Initialize pseudopotentials and local potentials
    #
    G2_shells = pw.gvec.G2_shells
    Ngl = length(G2_shells)
    Vgl = zeros(Float64, Ngl)
    idx_g2shells = pw.gvec.idx_g2shells
    V_Ps_loc = zeros(Float64, Npoints)
    Vg = zeros(ComplexF64, Npoints)


    #
    # Initialize pseudopotentials and local ionic potentials
    #
    is_psp_using_nlcc = zeros(Bool,Nspecies) # by default we don't use any NLCC    
    V_of_0 = 0.0
    #
    for isp in 1:Nspecies
        psp = pspots[isp] # shortcut
        if psp.is_nlcc
            is_psp_using_nlcc[isp] = true
        end
        eval_Vloc_G!( psp, G2_shells, Vgl )
        for ig in 1:Ng
            ip = idx_g2r[ig]
            igl = idx_g2shells[ig]
            Vg[ip] += strf[ig,isp] * Vgl[igl] / CellVolume
        end
    end
    V_of_0 += real(Vg[1]) # using PWSCF convention
    V_Ps_loc = real(G_to_R(pw, Vg))
    V_Ps_loc *= Npoints # Rescale using PWSCF convention


    if options.lspinorb
        Nspin_channel = 1
        Nspin_comp = 4
        @assert options.noncollinear
    end
    if options.noncollinear
        Nspin_channel = 1
        Nspin_comp = 4
    else
        Nspin_channel = options.Nspin_channel
        Nspin_comp = options.Nspin_comp
    end
    @info "Nspin_channel=$(Nspin_channel) Nspin_comp=$(Nspin_comp)"

    # other potential terms are set to zero
    V_Hartree = zeros( Float64, Npoints )
    V_xc = zeros( Float64, Npoints, Nspin_comp )
    V_loc_tot = zeros( Float64, Npoints, Nspin_comp )
    if pw.using_dual_grid
        # We initialize smooth local potential here (total)
        potentials = Potentials(
            V_Ps_loc, V_Hartree, V_xc, V_loc_tot,
            zeros(Float64, prod(pw.Nss), Nspin_comp),
            zeros(Float64, Npoints, Nspin_comp)
        )
    else
        potentials = Potentials(
            V_Ps_loc, V_Hartree, V_xc, V_loc_tot,
            nothing,
            zeros(Float64, Npoints, Nspin_comp)
        )
    end


    #
    energies = Energies()
    #
    rhoe = zeros( Float64, Npoints, Nspin_comp )
    #
    # Initialize core electron density for NLCC if needed
    #
    if any(is_psp_using_nlcc)
        rhoe_core = zeros(Float64, Npoints, 1) # FIXME: Nspin=1
        calc_rhoe_core!(atoms, pw, pspots, rhoe_core)
    else
        rhoe_core = nothing
    end

    # Initialize electronic states variable
    #
    # extra_states is given
    if options.extra_states > -1
        @info "Pass here 154"
        electrons = Electrons( atoms, pspots,
            Nspin_channel = Nspin_channel,
            Nkpt = kpoints.Nkpt,
            Nstates_empty = options.extra_states,
            noncollinear = options.noncollinear,
            domag = domag
        )
    # no extra_states is given but Nstates is given
    elseif options.Nstates > -1
        @info "Pass here 163"
        electrons = Electrons( atoms, pspots,
            Nspin_channel = Nspin_channel,
            Nkpt = kpoints.Nkpt,
            Nstates = options.Nstates,
            noncollinear = options.noncollinear,
            domag = domag
        )
    #
    elseif (options.Nstates == -1) && (options.extra_states == -1)
        @info "Pass here 172"
        # Default value for Nstates and Nstates_empty
        # Nstates will be calculated automatically
        electrons = Electrons( atoms, pspots,
            Nspin_channel = Nspin_channel,
            Nkpt = kpoints.Nkpt,
            noncollinear = options.noncollinear,
            domag = domag
        )
    else
        error("Error in initializing instance of Electrons")
    end


    # FIXME: Make parametric PsPotNL ?
    # NL pseudopotentials
    are_using_upf = zeros(Bool, Nspecies)
    for isp in 1:Nspecies
        pspfile = pspots[isp].pspfile
        are_using_upf[isp] = PWDFT.is_using_extension_upf(pspfile)
    end
    if all(are_using_upf)
        is_gga = (options.xcfunc == "PBE") # XXX FIX THIS !!!!
        pspotNL = PsPotNL_UPF( atoms, pw, pspots,
            is_gga = is_gga,
            Nspin = Nspin_comp,
            noncollinear = options.noncollinear,
            lspinorb = options.lspinorb
        )
    elseif all(.!are_using_upf)
        pspotNL = PsPotNL( atoms, pw, pspots, check_norm = false )
    else
        error("Not supporting mixed pseudopotential types: GTH and UPF")
    end
    # XXX For the moment we treat GTH and UPF as different.
    # We also do not support mixing UPF and GTH pspots
    # Ideally we should have one PsPotNL type only.


    atoms.Zvals = get_Zvals( pspots )

    ik = 1
    ispin = 1

    if sym_info.Nsyms > 1
        rhoe_symmetrizer = RhoeSymmetrizer( atoms, pw, sym_info )
    else
        rhoe_symmetrizer = RhoeSymmetrizer() # dummy rhoe_symmetrizer
    end

    # FIXME xc_calc constructor is not yet able to determine or check whether
    # given x_id or x_id is compatible with default or given is_gga or is_metagga
    # We force it here.
    if options.use_xc_internal
        if options.xcfunc == "SCAN" # This is broken
            xc_calc = XCCalculator(is_metagga=true) # id?
        elseif options.xcfunc == "PBE"
            xc_calc = XCCalculator(is_gga=true, x_id=101, c_id=130) # PBE
        else
            xc_calc = XCCalculator() # will use VWN
        end
    else
        # Using Libxc is the default
        if options.xcfunc == "SCAN"
            xc_calc = LibxcXCCalculator(is_metagga=true, Npoints=Npoints, Nspin=options.Nspin)
            # x_id and c_id are set to be 263 and 267 (SCAN), respectively
        elseif options.xcfunc == "PBE"
            xc_calc = LibxcXCCalculator(is_gga=true, x_id=101, c_id=130)
        else
            xc_calc = LibxcXCCalculator()
        end
    end


    need_overlap = any(pspotNL.are_ultrasoft) || any(pspotNL.are_paw)

    return Hamiltonian(
        pw, potentials, energies, rhoe, rhoe_core,
        electrons, atoms, sym_info, rhoe_symmetrizer,
        pspots, pspotNL, options.xcfunc, xc_calc, ik, ispin, need_overlap,
        options
    )
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
function Hamiltonian(
    atoms::Atoms,
    pspfiles::Vector{String},
    ecutwfc::Float64;
    # Keyword arguments
    dual::Float64=4.0,
    Nspin_channel::Int64=1,
    Nspin_comp::Int64=1,
    meshk::Vector{Int64}=[1,1,1],  # FIXME: convert to tuple?
    shiftk::Vector{Int64}=[0,0,0],
    time_reversal::Bool=true,
    Ns_::Tuple{Int64,Int64,Int64}=(0,0,0),
    kpoints::Union{KPoints,Nothing}=nothing,
    kpts_str::String="",
    xcfunc::String="VWN",
    use_xc_internal::Bool=false,
    extra_states::Int64=-1,
    Nstates::Int64=-1,
    use_symmetry::Bool=true,
    lspinorb::Bool=false,
    noncollinear::Bool=false,
    starting_magn = nothing,
    angle1 = nothing,
    angle2 = nothing,
    use_smearing = false
)
    # FIXME: This constructor is type-unstable (?)

    # Build HamiltonianOptions from kwargs
    options = HamiltonianOptions(
        dual, Nspin_channel, Nspin_comp, meshk, shiftk, time_reversal, Ns_,
        kpoints, kpts_str, xcfunc, use_xc_internal,
        extra_states, Nstates, use_symmetry, use_smearing,
        starting_magn, angle1, angle2,
        lspinorb, noncollinear
    )

    Nspecies = atoms.Nspecies
    if Nspecies != size(pspfiles)[1]
        error( @sprintf("Length of pspfiles is not equal to %d\n", Nspecies) )
    end

    #
    # XXX We don't support mixed pseudopotentials set
    #
    are_using_upf = zeros(Bool, Nspecies)
    for isp in 1:Nspecies
        are_using_upf[isp] = is_using_extension_upf(pspfiles[isp])
    end
    #
    if all(are_using_upf)
        pspots = Vector{PsPot_UPF}(undef,Nspecies)
        for isp in 1:Nspecies
            pspots[isp] = PsPot_UPF(pspfiles[isp])
        end
    elseif all(.!are_using_upf)
        pspots = Vector{PsPot_GTH}(undef,Nspecies)
        for isp in 1:Nspecies
            pspots[isp] = PsPot_GTH(pspfiles[isp])
        end
    else
        error("Not supporting mixed pseudopotential types: GTH and UPF")
    end

    return Hamiltonian(atoms, pspots, ecutwfc, options)
end


#=
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
    pspots = Array{PsPot_GTH}(undef,1)
    pspots[1] = PsPot_GTH()

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
    if pw.using_dual_grid
        # We initialize smooth local potential here (total)
        potentials = Potentials(
            V_Ps_loc, V_Hartree, V_xc, V_loc_tot,
            zeros(Float64, prod(pw.Nss), Nspin)
        )
    else
        potentials = Potentials(
            V_Ps_loc, V_Hartree, V_xc, V_loc_tot,
            nothing
        )
    end
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

    return Hamiltonian( pw, potentials, energies, rhoe, nothing,
                        electrons, atoms, sym_info, rhoe_symmetrizer,
                        pspots, pspotNL, xcfunc, xc_calc, ik, ispin )
end
=#


# FIXME: should be merged with Array{Float64,2} version.
# psiks is needed for metagga
function update!( Ham::Hamiltonian, psiks::BlochWavefunc, Rhoe::Vector{Float64} )

    # assumption Nspin = 1
    # Copy
    Ham.rhoe[:,1] .= Rhoe[:,1]
    
    pw = Ham.pw
    xc_calc = Ham.xc_calc
    V_Hartree = Ham.potentials.Hartree
    V_xc = Ham.potentials.XC
    V_Total = Ham.potentials.Total
    V_Ps_loc = Ham.potentials.Ps_loc

    Poisson_solve!(pw, Rhoe, V_Hartree)

    if !isnothing(Ham.rhoe_core)
        Rhoe[:] .+= Ham.rhoe_core
    end
    if Ham.xcfunc == "SCAN"
        @views calc_Vxc_SCAN!( Ham, psiks, Rhoe, V_xc[:,1] )
    #
    elseif Ham.xcfunc == "PBE"
        @views calc_Vxc_PBE!( xc_calc, pw, Rhoe, V_xc[:,1] )
    #
    else
        @views calc_Vxc_VWN!( xc_calc, Rhoe, V_xc[:,1] )
    end
    # Restore Rhoe
    if !isnothing(Ham.rhoe_core)
        Rhoe[:] .-= Ham.rhoe_core
    end
    #
    Npoints = size(Rhoe, 1)
    for ip in 1:Npoints
        V_Total[ip,1] = V_Ps_loc[ip] + V_Hartree[ip] + V_xc[ip,1]  
    end
    return
end




#
# TODO: save old total potential XXXXX
# XXX: Add option to save old potential to the Hamiltonian
#

"""
    update!(Ham, psiks, rhoe)

Update Ham.rhoe and calculate Hartree and XC potentials for given `rhoe` in real space.
"""
function update!(
    Ham::Hamiltonian,
    psiks::BlochWavefunc,
    Rhoe::Array{Float64,2}
)
    #
    Nspin = size(Rhoe, 2)
    Npoints = size(Rhoe, 1)
    if Nspin == 1
        update!(Ham, psiks, Rhoe[:,1])
        # No need for views, Rhoe is not updated
        return
    end
    # From this point on Nspin == 2
    @assert Nspin == 2

    # Some shortcuts
    V_H = Ham.potentials.Hartree
    V_xc = Ham.potentials.XC
    V_Ps_loc = Ham.potentials.Ps_loc
    #
    # Copy
    Ham.rhoe[:,:] .= Rhoe[:,:] # XXX need this? Or set this outside
    #
    Rhoe_total = Rhoe[:,1] + Rhoe[:,2] # Nspin is 2
    V_H[:] .= real( G_to_R( Ham.pw, Poisson_solve(Ham.pw, Rhoe_total) ) )
    # Use this instead?
    # Poisson_solve!(pw, Rhoe, V_Hartree)
    #
    # Add core density
    if !isnothing(Ham.rhoe_core)
        # Nspin == 2
        Rhoe[:,1] .+= Ham.rhoe_core*0.5
        Rhoe[:,2] .+= Ham.rhoe_core*0.5
    end
    #
    # FIXME: spinpol MetaGGA is not yet implemented
    if Ham.xcfunc == "PBE"
        calc_Vxc_PBE!( Ham.xc_calc, Ham.pw, Rhoe, V_xc )
    else 
        # VWN is the default
        calc_Vxc_VWN!( Ham.xc_calc, Rhoe, V_xc )
    end
    # Restore, because we don't modify the argument Rhoe
    if !isnothing(Ham.rhoe_core)
        # Nspin == 2
        Rhoe[:,1] .-= Ham.rhoe_core*0.5
        Rhoe[:,2] .-= Ham.rhoe_core*0.5
    end
    #
    # Set total potential
    for ispin in 1:Nspin, ip in 1:Npoints
        Ham.potentials.Total[ip,ispin] = V_Ps_loc[ip] + V_H[ip] + V_xc[ip,ispin]
    end
    return
end

# For compatibility
function update!( Ham::Hamiltonian, rhoe::Vector{Float64} )
    # assumption Nspin = 1
    Ham.rhoe[:,1] = rhoe
    Ham.potentials.Hartree = real( G_to_R( Ham.pw, Poisson_solve(Ham.pw, rhoe) ) )    
    if Ham.xcfunc == "PBE"
        Ham.potentials.XC[:,1] = calc_Vxc_PBE( Ham.xc_calc, Ham.pw, rhoe )
    else
        # VWN is the default
        if isnothing(Ham.rhoe_core)
            Ham.potentials.XC[:,1] = calc_Vxc_VWN( Ham.xc_calc, rhoe )
        else
            Ham.potentials.XC[:,1] = calc_Vxc_VWN( Ham.xc_calc, rhoe + Ham.rhoe_core )
        end
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