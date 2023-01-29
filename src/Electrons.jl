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
    Nspin::Int64
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
    Nspin = 1
    return Electrons( Nelectrons, Nstates, Nstates_occ, Focc, ebands, Nspin )
end

"""
Creates an instance of `Electrons` given following inputs:

- `atoms`: an instance of `Atoms`

- `Pspots`: an array of `PsPot_GTH` or `PsPot_UPF` instance(s)

- `Nspin`: (optional) number of spin. `Nspin=1` means without spin polarization.
  `Nspin=2` means with spin-polarization.

- `Nkpt`: (optional) number of kpoints

- `Nstates`: (optional) total number of electronic states

- `Nstates_empty`: (optional) number of additional states which will be
  regarded as empty.
"""
function Electrons(
    atoms::Atoms, pspots::Vector{T};
    Nspin=1, Nkpt=1,
    Nstates=-1, Nstates_empty=-1
) where T <: AbstractPsPot
    #
    return Electrons(
        atoms, get_Zvals(pspots), 
        Nspin=Nspin, Nkpt=Nkpt, Nstates=Nstates, Nstates_empty=Nstates_empty
    )
end


"""
Creates an instance of `Electrons` given following inputs:

- `atoms`: an instance of `Atoms`

- `zvals`: an array of `Float64` specifying number of valence electrons
  for each atomic species.

- `Nspin`: (optional) number of spin. `Nspin=1` means without spin polarization.
  `Nspin=2` means with spin-polarization.

- `Nkpt`: (optional) number of kpoints

- `Nstates`: (optional) total number of electronic states

- `Nstates_empty`: (optional) number of additional states which will be
  regarded as empty.
"""
function Electrons(
    atoms::Atoms, zvals::Vector{Float64};
    Nspin=1, Nkpt=1,
    Nstates=-1, Nstates_empty=-1
)

    @assert Nspin <= 2
    @assert length(zvals) == atoms.Nspecies

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
        end
    else
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
        end 
    end

    Focc = zeros(Float64,Nstates,Nkpt*Nspin)
    ebands = zeros(Float64,Nstates,Nkpt*Nspin)

    Nstates_occ = Nstates - Nstates_empty
    
    if Nspin == 1
        for ist = 1:Nstates_occ-1
            Focc[ist,:] .= 2.0
        end
        if is_odd
            Focc[Nstates_occ,:] .= 1.0
        else
            Focc[Nstates_occ,:] .= 2.0
        end
    else
        for ist = 1:Nstates_occ-1
            Focc[ist,:] .= 1.0
        end
        idx_up = 1:Nkpt
        if is_odd
            # assign only to the spin up
            Focc[Nstates_occ,idx_up] .= 1.0
        else
            Focc[Nstates_occ,:] .= 1.0
        end
    end

    _check_Focc(Focc, Nkpt, Nelectrons)

    return Electrons( Nelectrons, Nstates, Nstates_occ, Focc, ebands, Nspin )
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
function Electrons(
    atoms::Atoms, Pspots,
    NelectronsSpin::Tuple{Int64,Int64};
    Nkpt=1, Nstates_extra=0
)

    # Nstates_extra is always empty

    Nspin = 2
    Nelectrons = get_Nelectrons(atoms,Pspots)
    @assert round(Int64,Nelectrons) == sum(NelectronsSpin)

    Nstates_occ = maximum(NelectronsSpin)
    Nstates = Nstates_occ + Nstates_extra

    Focc = zeros(Float64,Nstates,Nkpt*Nspin)
    ebands = zeros(Float64,Nstates,Nkpt*Nspin)

    for ik in 1:Nkpt
        for i in 1:NelectronsSpin[1]
            Focc[i,ik] = 1.0
        end
        for i in 1:NelectronsSpin[2]
            Focc[i,Nkpt+ik] = 1.0
        end
    end

    _check_Focc(Focc, Nkpt, Nelectrons)

    return Electrons( Nelectrons, Nstates, Nstates_occ, Focc, ebands, Nspin )
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
