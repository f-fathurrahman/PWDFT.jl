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

include("Atoms_io.jl")
include("Atoms_utils.jl")

