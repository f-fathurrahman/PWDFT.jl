mutable struct Atoms
    Natoms::Int64
    Nspecies::Int64
    positions::Array{Float64,2}
    atm2species::Array{Int64,1}
    atsymbs::Array{String,1}
    SpeciesSymbols::Array{String,1}  # unique symbols
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
    @printf("Coordinates in bohr:\n")
    @printf("\n")
    for ia = 1:a.Natoms
        @printf("%5s %18.10f %18.10f %18.10f\n", a.atsymbs[ia],
                a.positions[1,ia], a.positions[2,ia], a.positions[3,ia])
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
    return Atoms( Natoms, Nspecies, positions, atm2species, atsymbs, SpeciesSymbols )
end


"""
Initialize from xyz file
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

    return Atoms(Natoms, Nspecies, positions, atm2species, atsymbs, SpeciesSymbols )

    #ALLOCATE( AtomicMasses(Nspecies) )
    #AtomicMasses(:) = 0.d0 ! FIXME: Use internal database to set this

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
