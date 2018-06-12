if VERSION > v"0.6.3"
    using Printf
    using LinearAlgebra
end

mutable struct Atoms
    Natoms::Int64
    Nspecies::Int64
    positions::Array{Float64,2}
    atm2species::Array{Int64,1}
    atsymbs::Array{String,1}         # for each atom
    SpeciesSymbols::Array{String,1}  # unique symbols
    LatVecs::Array{Float64,2}
    Zvals::Array{Float64,1}   # unique
end

# Overload println for Atoms
import Base.println
function println( a::Atoms )
    @printf("\n")
    @printf("                                     -----\n")
    @printf("                                     Atoms\n")
    @printf("                                     -----\n")
    @printf("\n")
    @printf("Natoms   = %5d\n", a.Natoms)
    @printf("Nspecies = %5d\n", a.Nspecies)
    @printf("\n")
    for isp = 1:a.Nspecies
        @printf("Species %3d: %s\n", isp, a.SpeciesSymbols[isp])
    end
    
    @printf("\n")
    @printf("Cartesian coordinates in bohr:\n")
    @printf("\n")
    for ia = 1:a.Natoms
        @printf("%5s %18.10f %18.10f %18.10f\n", a.atsymbs[ia],
                a.positions[1,ia], a.positions[2,ia], a.positions[3,ia])
    end
    
    @printf("\n")
    @printf("Direct lattice vectors: (v1, v2, and v3 are given by column)\n")
    @printf("\n")
    for i = 1:3
        @printf("%18.10f %18.10f %18.10f\n", a.LatVecs[i,1], a.LatVecs[i,2], a.LatVecs[i,3])
    end

    frac_pos = inv(a.LatVecs)*a.positions
    @printf("\n")
    @printf("Fractional coordinates:\n")
    @printf("\n")
    for ia = 1:a.Natoms
        @printf("%5s %18.10f %18.10f %18.10f\n", a.atsymbs[ia],
                frac_pos[1,ia], frac_pos[2,ia], frac_pos[3,ia])
    end

end

# dummy atoms, contains only one atom
function Atoms()
    Natoms = 1
    Nspecies = 1
    positions = zeros(3,Natoms)
    atm2species = [1]
    atsymbs = ["X"]
    SpeciesSymbols = ["X"]  # unique symbols
    LatVecs = 10.0*[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    Zvals = zeros(Nspecies)
    return Atoms( Natoms, Nspecies, positions, atm2species, atsymbs, SpeciesSymbols, LatVecs, Zvals )
end



function init_atoms_xyz_ext( filexyz; in_bohr=false, verbose=false )

end


"""
Initialize from a simple xyz file
"""
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
	#Nspecies = get_Nspecies( atsymbs )

    # Determine unique species symbols
	#SpeciesSymbols = get_SpeciesSymbols( Nspecies, atsymbs )

    SpeciesSymbols = unique(atsymbs)
    Nspecies = length(SpeciesSymbols)

    # Mapping of atoms to species index
    atm2species = get_atm2species( atsymbs, SpeciesSymbols )

    LatVecs = 10.0*[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    Zvals = zeros(Nspecies)
    return Atoms(Natoms, Nspecies, positions, atm2species, atsymbs, SpeciesSymbols, LatVecs, Zvals )

    #ALLOCATE( AtomicMasses(Nspecies) )
    #AtomicMasses(:) = 0.d0 ! FIXME: Use internal database to set this

end

"""
Just like `init_atoms_xyz` but instead of file, the contents are directly feed
to the function.
"""
function init_atoms_xyz_string(str::String; in_bohr=false, verbose=false)
    lines = split(str,"\n")
    l = lines[1]
    Natoms = parse(Int64, l)
    positions = zeros(3,Natoms)
    atsymbs = Array{String}(Natoms)
    #
    for ia = 1:Natoms
        ll = split( lines[2+ia] )
        atsymbs[ia] = ll[1]
        positions[1,ia] = parse( Float64, ll[2] )
        positions[2,ia] = parse( Float64, ll[3] )
        positions[3,ia] = parse( Float64, ll[4] )
    end

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
	#Nspecies = get_Nspecies( atsymbs )

    # Determine unique species symbols
	#SpeciesSymbols = get_SpeciesSymbols( Nspecies, atsymbs )

    SpeciesSymbols = unique(atsymbs)
    Nspecies = length(SpeciesSymbols)

    # Mapping of atoms to species index
    atm2species = get_atm2species( atsymbs, SpeciesSymbols )

    LatVecs = zeros(3,3)
    Zvals = zeros(Nspecies)
    return Atoms(Natoms, Nspecies, positions, atm2species, atsymbs, SpeciesSymbols, LatVecs, Zvals )
    
end

function get_Zatoms( atoms::Atoms )

    ZATOMS = Dict(
    "H"  => 1,
    "He" => 2,
    "Li" => 3,
    "Be" => 4,
    "B"  => 5,
    "C"  => 6,
    "N"  => 7,
    "O"  => 8,
    "F"  => 9,
    "Ne" => 10,
    "Na" => 11,
    "Mg" => 12,
    "Al" => 13,
    "Si" => 14,
    "P"  => 15,
    "S"  => 16,
    "Cl" => 17,
    "Ar" => 18,
    "K"  => 19,
    "Ca" => 20,
    "Sc" => 21,
    "Ti" => 22,
    "V"  => 23,
    "Cr" => 24,
    "Mn" => 25,
    "Fe" => 26,
    "Co" => 27,
    "Ni" => 28,
    "Cu" => 29,
    "Zn" => 30,
    "Ga" => 31,
    "Ge" => 32,
    "As" => 33,
    "Se" => 34,
    "Br" => 34,
    "Kr" => 36, )


    Nspecies = atoms.Nspecies
    SpeciesSymbols = atoms.SpeciesSymbols
    Zatoms = zeros(Float64,Nspecies)
    for isp = 1:Nspecies
        Zatoms[isp] = ZATOMS[SpeciesSymbols[isp]]
    end
    return Zatoms
end


##
## Helper functions
##

"""
Given all atomic symbols present in a system, `atsymbs`,
this function returns number of species present by searching for
unique symbols. This function use a naive algorithm.
Similar functionality can be obtained by using the built-in `unique` function.
```
Nspecies = length(unique(str))
```
"""
function get_Nspecies( atsymbs::Array{String,1} )
	# Determine number of species
    Nspecies = 0
	Natoms = size(atsymbs)[1]
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
	return Nspecies
end


"""
Given all atomic symbols present in a system, `atsymbs`,
this function returns unique atomic species symbol.
This function use a naive algorithm.
Similar functionality can be obtained by using the built-in `unique` function.
```
SpeciesSymbols = unique(str)
```
"""
function get_SpeciesSymbols( Nspecies, atsymbs )
    
	SpeciesSymbols = Array{String}(Nspecies)	
	Natoms = size(atsymbs)[1]

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

    return SpeciesSymbols

end

"""
Get `atm2species` which described mapping between `atsymbs` and `SpeciesSymbols`.
"""
function get_atm2species( atsymbs, SpeciesSymbols )

    Natoms = size(atsymbs)[1]
    Nspecies = size(SpeciesSymbols)[1]

    atm2species = Array{Int64}(Natoms)

    for ia = 1:Natoms
        for isp = 1:Nspecies
            if atsymbs[ia] == SpeciesSymbols[isp]
                atm2species[ia] = isp
            end
        end
    end

    return atm2species
end
