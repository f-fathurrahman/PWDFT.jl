mutable struct Atoms
    Natoms::Int64
    Nspecies::Int64
    positions::Array{Float64,2}
    atm2species::Array{Int64,1}
    atsymbs::Array{String,1}
    SpeciesSymbols::Array{String,1}  # unique symbols
    Zvals::Array{Float64,1}
end

# Overload println for Atoms
import Base.println
function println( a::Atoms )
    @printf("\nAtoms information\n\n")
    @printf("Natoms   = %5d\n", a.Natoms)
    @printf("Nspecies = %5d\n", a.Nspecies)
    @printf("\nCoordinates in bohr:\n\n")
    for ia = 1:a.Natoms
        @printf("%5s %18.10f %18.10f %18.10f\n", a.atsymbs[ia],
                a.positions[1,ia], a.positions[2,ia], a.positions[3,ia])
    end
end


function init_atoms_xyz(filexyz; in_bohr=false, verbose=false)
    f = open(filexyz, "r")
    l = readline(f)
    Natoms = parse(Int64, l)
    positions = zeros(3,Natoms)
    atsymbs = Array{String}(Natoms)
    l = readline(f)
    for ia = 1:Natoms
        ll = split( readline(f) )
        atsymbs[ia] = ll[1]
        positions[1,ia] = parse( Float64, ll[2] )
        positions[2,ia] = parse( Float64, ll[3] )
        positions[3,ia] = parse( Float64, ll[4] )
    end
    close(f)

    # convert from angstrom to bohr
    if !in_bohr
        if verbose
            println("Coordinate in xyz file is assumed to be given in angstrom")
            println("It will be converted to bohr")
        end
        positions[:,:] = positions[:,:] * ANG2BOHR
    else
        if verbose
            println("Coordinate in xyz file is assumed to be given in bohr")
        end
    end

    # Determine number of species
    Nspecies = 0
    for ia = 1:Natoms
        k2 = 0
        for k1 = 1:ia-1
            if atsymbs[k1] == atsymbs[ia]
                k2 = 1
            end
        end
        # find different
        if k2 == 0
            Nspecies = Nspecies + 1
        end
    end

    SpeciesSymbols = Array{String}(Nspecies)

    idx1 = 0
    for ia = 1:Natoms
        k2 = 0
        for k1 = 1:ia-1
            if atsymbs[k1] == atsymbs[ia]
                k2 = 1
            end
        end
        # Found different species
        if k2==0
            idx1 = idx1 + 1
            SpeciesSymbols[idx1] = atsymbs[ia]
        end
    end

    # Mapping of atoms to species index
    atm2species = Array{Int64}(Natoms)
    for ia = 1:Natoms
        for isp = 1:Nspecies
            if atsymbs[ia] == SpeciesSymbols[isp]
                atm2species[ia] = isp
            end
        end
    end

    # Set to zeros
    # It will be set later if using
    Zvals = zeros(Float64,Nspecies)

    return Atoms(Natoms, Nspecies, positions, atm2species, atsymbs, SpeciesSymbols, Zvals )

    #ALLOCATE( AtomicMasses(Nspecies) )
    #AtomicMasses(:) = 0.d0 ! FIXME: Use internal database to set this

end
