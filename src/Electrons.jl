"""
The type for describing electronic variables such as number of
electrons, occupation numbers, and energy levels.
"""
mutable struct Electrons
    Nelectrons::Float64
    Nstates::Int64
    Nstates_occ::Int64
    Focc::Array{Float64,2}
    ebands::Array{Float64,2}
    Nspin_channel::Int64
    use_smearing::Bool
    kT::Float64
    noncolin::Bool
    E_fermi::Float64
    Nspin_comp::Int64
end

"""
Creates a 'dummy' instance of `Electrons` with only one electron.
"""
function Electrons()
    Nelectrons = 1
    Nstates = 1
    Nstates_occ = 1
    Focc = zeros(Nstates,1) # Nkpt=1
    ebands = zeros(Nstates,1) # use Nkpt=1
    Nspin_channel = 1
    Nspin_comp = 1
    use_smearing = false
    kT = 0.0
    noncolin = false
    E_fermi = 0.0
    return Electrons(
        Nelectrons, Nstates, Nstates_occ, Focc, ebands, Nspin_channel,
        use_smearing, kT, noncolin, E_fermi, Nspin_comp
    )
end

"""
Creates an instance of `Electrons` given following inputs:

- `atoms`: an instance of `Atoms`

- `Pspots`: an array of `PsPot_GTH` or `PsPot_UPF` instance(s)

- `Nspin_channel`: (optional) number of spin. `Nspin_channel=1` means without spin polarization.
  `Nspin_channel=2` means with spin-polarization.

- `Nkpt`: (optional) number of kpoints

- `Nstates`: (optional) total number of electronic states

- `Nstates_empty`: (optional) number of additional states which will be
  regarded as empty.
"""
function Electrons(
    atoms::Atoms, pspots::Vector{T};
    Nspin_channel=1, Nkpt=1,
    Nstates=-1, Nstates_empty=-1, noncolin=false
) where T <: AbstractPsPot
    #
    return Electrons(
        atoms, get_Zvals(pspots), 
        Nspin_channel=Nspin_channel, Nkpt=Nkpt, Nstates=Nstates, Nstates_empty=Nstates_empty,
        noncolin=noncolin
    )
end


"""
Creates an instance of `Electrons` given following inputs:

- `atoms`: an instance of `Atoms`

- `zvals`: an array of `Float64` specifying number of valence electrons
  for each atomic species.

- `Nspin_channel`: (optional) number of spin. `Nspin_channel=1` means without spin polarization.
  `Nspin_channel=2` means with spin-polarization.

- `Nkpt`: (optional) number of kpoints

- `Nstates`: (optional) total number of electronic states

- `Nstates_empty`: (optional) number of additional states which will be
  regarded as empty.
"""
function Electrons(
    atoms::Atoms, zvals::Vector{Float64};
    Nspin_channel=1, Nkpt=1,
    Nstates=-1, Nstates_empty=-1, noncolin=false
)
    if !noncolin
        @assert Nspin_channel <= 2
    end
    @assert length(zvals) == atoms.Nspecies

    # Determine Nspin_comp
    Nspin_comp = 1 # default
    # Collinear magnetism
    if !noncolin && (Nspin_channel == 2)
        Nspin_comp = 2
    end
    # For noncolin Nspin_channel = 1 and Nspin_comp = 4
    if noncolin
        @assert Nspin_channel == 1
        Nspin_comp = 4
    end

    Nelectrons = get_Nelectrons(atoms, zvals)

    is_odd = round(Int64,Nelectrons)%2 == 1

    # If Nstates is not specified and Nstates_empty == 0, we calculate
    # Nstates manually from Nelectrons
    if Nstates == -1
        Nstates = round(Int64, Nelectrons/2)
        if Nstates*2 < Nelectrons
            Nstates = Nstates + 1
        end
        if Nstates_empty > 0
            Nstates = Nstates + Nstates_empty
        else
            Nstates_empty = 0 # Given Nstates_empty a valid value
        end
    else
        # Nstates is not given its default (invalid value)
        if Nstates_empty != -1
            println("Please specify Nstates only or Nstates_empty only")
            error()
        end
        NstatesMin = round(Int64, Nelectrons/2)
        if NstatesMin*2 < Nelectrons
            NstatesMin += 1
        end
        if Nstates < NstatesMin
           @printf("Given Nstates is not enough: Nstates = %d, NstatesMin = %d", Nstates, NstatesMin)
        end
        if Nstates > NstatesMin
            Nstates_empty = Nstates - NstatesMin
        else
            Nstates_empty = 0
        end
    end
    # Nstates, Nstates_empty must be set up to their valid values now
    # i.e. they should not have value of -1

    #println("Nstates = ", Nstates)
    #println("Nstates_empty = ", Nstates_empty)

    Focc = zeros(Float64,Nstates,Nkpt*Nspin_channel)
    ebands = zeros(Float64,Nstates,Nkpt*Nspin_channel)
    
    Nstates_occ = Nstates - Nstates_empty

    if noncolin
        OCC_MAX = 1.0
    else
        if Nspin_channel == 1
            OCC_MAX = 2.0
        else
            OCC_MAX = 1.0
        end
    end

    if Nspin_channel == 1
        for ist in 1:Nstates_occ-1
            Focc[ist,:] .= OCC_MAX
        end
        if is_odd
            Focc[Nstates_occ,:] .= 1.0
        else
            Focc[Nstates_occ,:] .= 2.0
        end
    else
        for ist = 1:Nstates_occ-1
            Focc[ist,:] .= OCC_MAX
        end
        idx_up = 1:Nkpt
        if is_odd
            # assign only to the spin up
            Focc[Nstates_occ,idx_up] .= OCC_MAX
        else
            Focc[Nstates_occ,:] .= OCC_MAX
        end
    end

    _check_Focc(Focc, Nkpt, Nelectrons)

    use_smearing = false
    kT = 0.0
    E_fermi = 0.0
    return Electrons(
        Nelectrons, Nstates, Nstates_occ, Focc, ebands, Nspin_channel,
        use_smearing, kT, noncolin, E_fermi, Nspin_comp
    )
end


function _check_Focc(Focc::Matrix{Float64}, Nkpt::Int64, Nelectrons)
    sFocc = sum(Focc)/Nkpt
    # Check if the generated Focc is consistent
    if abs( sFocc - Nelectrons ) > eps()
        @printf("sum Focc = %f, Nelectrons = %f\n", sFocc, Nelectrons)
        error(@sprintf("ERROR diff sum(Focc) and Nelectrons is not small\n"))
    end
    return
end


"""
NelectronsSpin = (Nel_up, Nel_dn)
This is useful for molecule.
"""
# XXX: Not tested for noncolin
function Electrons(
    atoms::Atoms, Pspots,
    NelectronsSpin::Tuple{Int64,Int64};
    Nkpt=1, Nstates_extra=0
)

    # Nstates_extra is always empty

    Nspin_channel = 2
    Nelectrons = get_Nelectrons(atoms,Pspots)
    @assert round(Int64,Nelectrons) == sum(NelectronsSpin)

    Nstates_occ = maximum(NelectronsSpin)
    Nstates = Nstates_occ + Nstates_extra

    Focc = zeros(Float64,Nstates,Nkpt*Nspin_channel)
    ebands = zeros(Float64,Nstates,Nkpt*Nspin_channel)

    for ik in 1:Nkpt
        for i in 1:NelectronsSpin[1]
            Focc[i,ik] = 1.0
        end
        for i in 1:NelectronsSpin[2]
            Focc[i,Nkpt+ik] = 1.0
        end
    end

    _check_Focc(Focc, Nkpt, Nelectrons)

    use_smearing = false
    kT = 0.0
    E_fermi = 0.0
    noncolin = false
    Nspin_comp = 2 # collinear magnetism
    return Electrons(
        Nelectrons, Nstates, Nstates_occ, Focc, ebands, Nspin_channel,
        use_smearing, kT, noncolin, E_fermi
    )
end


"""
Returns number of electrons for a given `atoms::Atoms` and
`Pspots::Array{AbstractPsPot,1}`. Number of electrons will be
calculated as sum of valence electrons for each atom.
"""
function get_Nelectrons( atoms::Atoms, pspots::Vector{T} ) where T <: AbstractPsPot
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    Nelectrons = 0.0
    for ia in 1:Natoms
        isp = atm2species[ia]
        Nelectrons += pspots[isp].zval
    end
    return Nelectrons
end



function get_Nelectrons( atoms::Atoms, zvals::Vector{Float64} ) 
    @assert length(zvals) == atoms.Nspecies
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    Nelectrons = 0.0
    for ia in 1:Natoms
        isp = atm2species[ia]
        Nelectrons += zvals[isp]
    end
    return Nelectrons
end



"""
Returns array `Zvals[1:Nspecies]` from a given `PsPots`.
"""
function get_Zvals( pspots::Vector{T} ) where T <: AbstractPsPot
    Nspecies = length(pspots)
    zvals = zeros(Float64, Nspecies)
    for isp in 1:Nspecies
        zvals[isp] = pspots[isp].zval
    end
    return zvals
end

include("Electrons_io.jl")
