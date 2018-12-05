"""
The type for representing an atomic-scale structures such as molecules
and solids.
"""
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


"""
Creates an instance of `Atoms`.

If no arguments are given, then a dummy `Atoms` is returned.

Use only one of `xyz_file`, `xyz_string` and `xyz_string_frac` for specifying
number of atoms, atomic symbols, and positions, following xyz format,
for example:
```
3

C   0.0   0.0   0.0
O   0.0   0.0   1.5
```

The argument `xyz_file` is used when the structure is given in an xyz file, which
in this case `xyz_file` is a string of path to the xyz file. By default the unit
used is angstrom. If the coordinate is given in bohr, set `in_bohr=true`.

Alternatively, one can pass the content of xyz file as string in `xyz_string`.

`xyz_string_frac` is the same as `xyz_string`, however the actual coordinates will be
transformed by multiplying it with `LatVecs`. This is useful for crystals.
"""
function Atoms( ;xyz_file="", xyz_string="", xyz_string_frac="", ext_xyz_file="",
    in_bohr=false, LatVecs=10*[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0] )

    if xyz_file != ""
        atoms = init_atoms_xyz(xyz_file, in_bohr=in_bohr)
        atoms.LatVecs = LatVecs
        return atoms
    elseif xyz_string != ""
        atoms = init_atoms_xyz_string(xyz_string, in_bohr=in_bohr)
        atoms.LatVecs = LatVecs
        return atoms
    elseif xyz_string_frac != ""
        atoms = init_atoms_xyz_string(xyz_string_frac, in_bohr=in_bohr)
        atoms.positions = LatVecs*atoms.positions
        atoms.LatVecs = LatVecs
        return atoms
    elseif ext_xyz_file != ""
        atoms = init_atoms_xyz_ext(ext_xyz_file, in_bohr=in_bohr)
        return atoms
    else
        # No arguments are assumed to be provided
        # dummy `atoms`, contains only one atom is returned
        Natoms = 1
        Nspecies = 1
        positions = zeros(3,Natoms)
        atm2species = [1]
        atsymbs = ["X"]
        SpeciesSymbols = ["X"]  # unique symbols
        Zvals = zeros(Nspecies)
        return Atoms( Natoms, Nspecies, positions, atm2species, atsymbs, SpeciesSymbols, LatVecs, Zvals )
    end

end


# extended XYZ format used in ASE
function init_atoms_xyz_ext( filexyz; in_bohr=false, verbose=false )
    f = open(filexyz, "r")
    l = readline(f)
    Natoms = parse(Int64, l)
    positions = zeros(3,Natoms)
    atsymbs = Array{String}(undef,Natoms)
    #
    l = readline(f)
    ll = replace( split(l,"=")[2], "\"" => "") # remove "s
    rll = split(ll)
    LatVecs_ = zeros(9)
    for i = 1:9
        LatVecs_[i] = parse( Float64, rll[i] )
    end
    LatVecs = reshape(LatVecs_, (3,3))
    #
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
            println("Coordinates and lattice vectors in extended xyz file is assumed")
            println("to be given in angstrom")
            println("It will be converted to bohr")
        end
        positions[:,:] = positions[:,:] * ANG2BOHR
        LatVecs[:,:] = LatVecs[:,:] * ANG2BOHR
    else
        if verbose
            println("Coordinate in xyz file is assumed to be given in bohr")
        end
    end

    SpeciesSymbols = unique(atsymbs)
    Nspecies = length(SpeciesSymbols)

    # Mapping of atoms to species index
    atm2species = get_atm2species( atsymbs, SpeciesSymbols )

    Zvals = zeros(Nspecies)
    return Atoms(Natoms, Nspecies, positions, atm2species, atsymbs, SpeciesSymbols, LatVecs, Zvals )

end


"""
Creates an instance of `Atoms` from a simple xyz file.
"""
function init_atoms_xyz(xyz_file; in_bohr=false, verbose=false)
    f = open(xyz_file, "r")
    l = readline(f)
    Natoms = parse(Int64, l)
    positions = zeros(3,Natoms)
    atsymbs = Array{String}(undef,Natoms)
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

    SpeciesSymbols = unique(atsymbs)
    Nspecies = length(SpeciesSymbols)

    # Mapping of atoms to species index
    atm2species = get_atm2species( atsymbs, SpeciesSymbols )

    LatVecs = 10.0*[1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    Zvals = zeros(Nspecies)
    return Atoms(Natoms, Nspecies, positions, atm2species, atsymbs, SpeciesSymbols, LatVecs, Zvals )

end

"""
Just like `init_atoms_xyz` but instead of file, the contents (of type String)
are directly feed to the function.
"""
function init_atoms_xyz_string(str::String; in_bohr=false, verbose=false)
    lines = split(str,"\n")
    l = lines[1]
    Natoms = parse(Int64, l)
    positions = zeros(3,Natoms)
    atsymbs = Array{String}(undef,Natoms)
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

    SpeciesSymbols = unique(atsymbs)
    Nspecies = length(SpeciesSymbols)

    # Mapping of atoms to species index
    atm2species = get_atm2species( atsymbs, SpeciesSymbols )

    LatVecs = zeros(3,3)
    Zvals = zeros(Nspecies)
    return Atoms(Natoms, Nspecies, positions, atm2species, atsymbs, SpeciesSymbols, LatVecs, Zvals )
    
end

# Overload println for Atoms
import Base: println

function println( a::Atoms; header=true )
    if header
        @printf("\n")
        @printf("                                     -----\n")
        @printf("                                     Atoms\n")
        @printf("                                     -----\n")
        @printf("\n")
    end
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

const ZATOMS = Dict(
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
    "Br" => 35,
    "Kr" => 36,
    "Rb" => 37,
    "Sr" => 38,
    "Y"  => 39,
    "Zr" => 40,
    "Nb" => 41,
    "Mo" => 42,
    "Tc" => 43,
    "Ru" => 44,
    "Rh" => 45,
    "Pd" => 46,
    "Ag" => 47,
    "Cd" => 48,
    "In" => 49,
    "Sn" => 50,
    "Sb" => 51,
    "Te" => 52,
    "I"  => 53,
    "Xe" => 54,
    "Cs" => 55,
    "Ba" => 56,
    "La" => 57,
    "Ce" => 58,
    "Pr" => 59,
    "Nd" => 60,
    "Pm" => 61,
    "Sm" => 62,
    "Eu" => 63,
    "Gd" => 64,
    "Tb" => 65,
    "Dy" => 66,
    "Ho" => 67,
    "Er" => 68,
    "Tm" => 69,
    "Yb" => 70,
    "Lu" => 71,
    "Hf" => 72,
    "Ta" => 73,
    "W"  => 74,
    "Re" => 75,
    "Os" => 76,
    "Ir" => 77,
    "Pt" => 78,
    "Au" => 79,
    "Hg" => 80,
    "Tl" => 81,
    "Pb" => 82,
    "Bi" => 83,
    "Po" => 84,
    "At" => 85,
    "Rn" => 86,)


"""
Returns an array of atomic numbers (with size `Nspecies`) given an instance
of `Atoms`.
"""
function get_Zatoms( atoms::Atoms )
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
    
	SpeciesSymbols = Array{String}(undef,Nspecies)
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

    atm2species = Array{Int64}(undef,Natoms)

    for ia = 1:Natoms
        for isp = 1:Nspecies
            if atsymbs[ia] == SpeciesSymbols[isp]
                atm2species[ia] = isp
            end
        end
    end

    return atm2species
end
