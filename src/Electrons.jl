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
function Electrons( atoms::Atoms, Pspots;
                    Nspin=1, Nkpt=1,
                    Nstates=nothing, Nstates_empty=0 )
    
    @assert( Nspin <= 2 )

    Nelectrons = get_Nelectrons(atoms,Pspots)

    is_odd = round(Int64,Nelectrons)%2 == 1

    # If Nstates is not specified and Nstates_empty == 0, we calculate
    # Nstates manually from Nelectrons
    if (Nstates == nothing)
        Nstates = round( Int64, Nelectrons/2 )
        if Nstates*2 < Nelectrons
            Nstates = Nstates + 1
        end
        if Nstates_empty > 0
            Nstates = Nstates + Nstates_empty
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

    sFocc = sum(Focc)/Nkpt
    # Check if the generated Focc is consistent
    if abs( sFocc - Nelectrons ) > eps()
        @printf("sum Focc = %f, Nelectrons = %f\n", sFocc, Nelectrons)
        error(@sprintf("ERROR: diff sum(Focc) and Nelectrons is not small\n"))
    end

    return Electrons( Nelectrons, Nstates, Nstates_occ, Focc, ebands, Nspin )
end

"""
Creates an instance of `Electrons` given following inputs:

- `atoms`: an instance of `Atoms`

- `Zvals`: an array of `Float64` specifying number of valence electrons
  for each atomic species.

- `Nspin`: (optional) number of spin. `Nspin=1` means without spin polarization.
  `Nspin=2` means with spin-polarization.

- `Nkpt`: (optional) number of kpoints

- `Nstates`: (optional) total number of electronic states

- `Nstates_empty`: (optional) number of additional states which will be
  regarded as empty.
"""
function Electrons( atoms::Atoms, Zvals::Array{Float64,1};
                    Nspin=1, Nkpt=1,
                    Nstates=nothing, Nstates_empty=0 )

    @assert( Nspin <= 2 )

    Nelectrons = 0.0
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    for ia = 1:Natoms
        isp = atm2species[ia]
        Nelectrons = Nelectrons + Zvals[isp]
    end

    is_odd = round(Int64,Nelectrons)%2 == 1

    # If Nstates is not specified and Nstates_empty == 0, we calculate
    # Nstates manually from Nelectrons
    if (Nstates == nothing)
        Nstates = round( Int64, Nelectrons/2 )
        if Nstates*2 < Nelectrons
            Nstates = Nstates + 1
        end
        if Nstates_empty > 0
            Nstates = Nstates + Nstates_empty
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

    sFocc = sum(Focc)/Nkpt
    # Check if the generated Focc is consistent
    if abs( sFocc - Nelectrons ) > eps()
        @printf("sum Focc = %f, Nelectrons = %f\n", sFocc, Nelectrons)
        error(@sprintf("ERROR diff sum(Focc) and Nelectrons is not small\n"))
    end

    return Electrons( Nelectrons, Nstates, Nstates_occ, Focc, ebands, Nspin )
end

"""
NelectronsSpin = (Nel_up, Nel_dn)
"""
function Electrons(
    atoms::Atoms, Pspots,
    NelectronsSpin::Tuple{Int64,Int64};
    Nkpt=1, Nstates_extra=0
)
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

    sFocc = sum(Focc)/Nkpt
    # Check if the generated Focc is consistent
    if abs( sFocc - Nelectrons ) > eps()
        @printf("sum Focc = %f, Nelectrons = %f\n", sFocc, Nelectrons)
        error(@sprintf("ERROR: diff sum(Focc) and Nelectrons is not small\n"))
    end

    return Electrons( Nelectrons, Nstates, Nstates_occ, Focc, ebands, Nspin )
end


"""
Returns number of electrons for a given `atoms::Atoms` and
`Pspots::Array{AbstractPsPot,1}`. Number of electrons will be
calculated as sum of valence electrons for each atom.
"""
function get_Nelectrons( atoms::Atoms, Pspots )
    Nelectrons = 0.0
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    for ia = 1:Natoms
        isp = atm2species[ia]
        Nelectrons = Nelectrons + Pspots[isp].zval
    end
    return Nelectrons
end

"""
Returns array `Zvals[1:Nspecies]` from a given `PsPots`.
"""
function get_Zvals( PsPots )
    Nspecies = size(PsPots)[1]
    Zvals = zeros(Float64, Nspecies)
    for isp = 1:Nspecies
        Zvals[isp] = PsPots[isp].zval
    end
    return Zvals
end

include("Electrons_io.jl")
