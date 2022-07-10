mutable struct Hamiltonian{Txc<:AbstractXCCalculator,Tpsp<:AbstractPsPot}
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
    pspotNL::Union{PsPotNL,PsPotNL_UPF}
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
function Hamiltonian(
    atoms::Atoms,
    pspfiles::Vector{String},
    ecutwfc::Float64;
    # Keyword arguments
    Nspin::Int64=1,
    meshk::Vector{Int64}=[1,1,1],
    shiftk::Vector{Int64}=[0,0,0],
    time_reversal::Bool=true,
    Ns_::Tuple{Int64,Int64,Int64}=(0,0,0),
    kpoints::Union{KPoints,Nothing}=nothing,
    kpts_str::String="",
    xcfunc::String="VWN",
    use_xc_internal::Bool=false,
    extra_states::Int64=0,
    use_symmetry::Bool=true
)

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
    G2_shells = pw.gvec.G2_shells
    Ngl = length(G2_shells)
    Vgl = zeros(Float64, Ngl)
    idx_g2shells = pw.gvec.idx_g2shells
    V_Ps_loc = zeros(Float64, Npoints)
    Vg = zeros(ComplexF64, Npoints)

    # XXX We don't support mixed pseudopotentials set
    if is_using_extension_upf(pspfiles[1])
        pspots = Array{PsPot_UPF}(undef,Nspecies)
    else
        pspots = Array{PsPot_GTH}(undef,Nspecies)
    end    

    #
    # Initialize pseudopotentials and local ionic potentials
    #
    is_psp_using_nlcc = zeros(Bool,Nspecies) # by default we don't use any NLCC
    #
    for isp = 1:Nspecies
        #
        if is_using_extension_upf(pspfiles[isp])            
            #
            # Using UPF
            #
            println("\nUsing UPF\n")
            pspots[isp] = PsPot_UPF( pspfiles[isp] )
            # build the interpolation table needed for nonlocal projectors
            _build_prj_interp_table!( pspots[isp], pw )
            #
            if pspots[isp].is_nlcc
                is_psp_using_nlcc[isp] = true
                println()
                @printf("!!! WARNING: %s is generated with nonlinear-core correction\n",
                    pspfiles[isp])
                @printf("!!! WARNING: This is not yet supported\n")
                println()
            end
        else
            #
            # Using GTH pseudopotential
            #
            pspots[isp] = PsPot_GTH( pspfiles[isp] )
        end
        #
        psp = pspots[isp] # shortcut
        #
        eval_Vloc_G!( psp, G2_shells, Vgl )
        for ig = 1:Ng
            ip = idx_g2r[ig]
            igl = idx_g2shells[ig]
            Vg[ip] = strf[ig,isp] * Vgl[igl]
        end
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
    #
    # Initialize core electron density for NLCC if needed
    #
    if any(is_psp_using_nlcc)
        rhoe_core = zeros(Float64, Npoints, 1) # FIXME: Nspin=1
        calc_rhoe_core!(atoms, pw, pspots, rhoe_core)
        println("sum rhoe_core = ", sum(rhoe_core))
    else
        rhoe_core = nothing
    end

    electrons = Electrons( atoms, pspots, Nspin=Nspin, Nkpt=kpoints.Nkpt,
                           Nstates_empty=extra_states )

    # NL pseudopotentials
    are_using_upf = zeros(Bool, Nspecies)
    for isp in 1:Nspecies
        are_using_upf[isp] = is_using_extension_upf(pspfiles[isp])
    end
    if all(are_using_upf)
        pspotNL = PsPotNL_UPF(atoms, pw, pspots)
    elseif all(.!are_using_upf)
        pspotNL = PsPotNL( atoms, pw, pspots, check_norm=false )
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

    return Hamiltonian( pw, potentials, energies, rhoe, rhoe_core,
                        electrons, atoms, sym_info, rhoe_symmetrizer,
                        pspots, pspotNL, xcfunc, xc_calc, ik, ispin )
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

    return Hamiltonian( pw, potentials, energies, rhoe, nothing,
                        electrons, atoms, sym_info, rhoe_symmetrizer,
                        pspots, pspotNL, xcfunc, xc_calc, ik, ispin )
end

# FIXME: should be merged with Array{Float64,2} version.
# psiks is needed for metagga
function update!( Ham::Hamiltonian, psiks::BlochWavefunc, rhoe::Array{Float64,1} )

    # assumption Nspin = 1
    Ham.rhoe[:,1] = rhoe
    
    pw = Ham.pw
    xc_calc = Ham.xc_calc
    V_Hartree = Ham.potentials.Hartree
    V_XC = Ham.potentials.XC
    V_Total = Ham.potentials.Total
    V_Ps_loc = Ham.potentials.Ps_loc

    Poisson_solve!(pw, rhoe, V_Hartree)

    if Ham.xcfunc == "SCAN"
        Vxc_tmp = zeros(size(rhoe,1)) # FIXME: use V_XC directly
        calc_Vxc_SCAN!( xc_calc, pw, psiks, rhoe, Vxc_tmp )
        V_XC[:,1] = Vxc_tmp[:]
    #
    elseif Ham.xcfunc == "PBE"
        V_XC[:,1] = calc_Vxc_PBE( xc_calc, pw, rhoe )
    #
    else
        # VWN is the default
        if Ham.rhoe_core == nothing
            @views V_XC[:,1] = calc_Vxc_VWN( xc_calc, rhoe )
        else
            @views V_XC[:,1] = calc_Vxc_VWN( xc_calc, rhoe + Ham.rhoe_core )
        end
    end
    
    Npoints = size(rhoe,1)
    for ip in 1:Npoints
        V_Total[ip,1] = V_Ps_loc[ip] + V_Hartree[ip] + V_XC[ip,1]  
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
        # FIXME: NLCC is not yet handled
    else 
        # VWN is the default
        if Ham.rhoe_core == nothing
            Ham.potentials.XC = calc_Vxc_VWN( Ham.xc_calc, rhoe )
        else
            Ham.potentials.XC = calc_Vxc_VWN( Ham.xc_calc, rhoe + Ham.rhoe_core )
        end
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
    else
        # VWN is the default
        if Ham.rhoe_core == nothing
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