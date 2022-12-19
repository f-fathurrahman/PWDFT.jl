struct PWSCFInput
    atoms::Atoms
    ecutwfc::Float64
    ecutrho::Float64
    pspfiles::Vector{String}
    meshk::Tuple{Int64,Int64,Int64}
end


#
# A naive function to read pwscf input
# Only a subset of possible combinations of parameters are considered.
#
function PWSCFInput( filename::String )

    # Default values, some are not valid for PWSCF
    acell = -1.0

    Natoms = 0
    Nspecies = 0

    ecutwfc = -1.0
    ecutrho = -1.0

    LatVecs = zeros(3,3)
    is_parse_cell = false
    N_parse_cell = 0

    xyz_string_frac = ""
    xyz_string_angstrom = ""
    xyz_string_bohr = ""
    is_parse_xyz = false
    N_parse_xyz = 0

    is_parse_kpoints = false
    N_parse_kpoints = 0
    meshk1 = 0
    meshk2 = 0
    meshk3 = 0

    in_angstrom = false
    in_fraction = false
    in_bohr = false

    is_parse_species = false
    N_parse_species = 0
    pspfiles = String[]
    species_symbols = String[]
    species_masses = Float64[]
    pseudo_dir = "./"

    f = open(filename, "r")
    
    while !eof(f)
        
        l = readline(f)
        
        # FIXME: This is not robust
        if occursin("  A =", l)
            ll = split(l, "=")
            acell = parse(Float64,ll[end])*ANG2BOHR
        end

        # Read number of atoms
        if occursin("nat =", l)
            ll = split(l, "=", keepempty=false)
            Natoms = parse(Int64, ll[end])
            println("Read Natoms = ", Natoms)
        end

        # Read number of species (ntyp)
        if occursin("ntyp =", l)
            ll = split(l, "=", keepempty=false)
            Nspecies = parse(Int64, ll[end])
            println("Read Nspecies = ", Nspecies)
        end


        # Read ecutwfc
        if occursin("ecutwfc =", l)
            ll = split(l, "=", keepempty=false)
            ecutwfc = parse(Float64, ll[end])
            println("Read ecutwfc (in Ry) = ", ecutwfc)
        end

        # Read ecutrho
        if occursin("ecutrho =", l)
            ll = split(l, "=", keepempty=false)
            ecutrho = parse(Float64, ll[end])
            println("Read ecutrho (in Ry) = ", ecutrho)
        end

        # Read pseudo_dir
        if occursin("pseudo_dir =", l)
            # FIXME: Notice the spaces around the equal sign
            # This is to prevent space before / in the path string
            ll = split(l, " = ", keepempty=false)
            pseudo_dir = replace(replace(ll[end], "'" => ""), "\"" => "")
            println("Read pseudo_dir = ", pseudo_dir)
        end


        #
        # Read cell parameters
        #
        # FIXME: They are assumed to be given in bohr !!!!
        if occursin("CELL_PARAMETERS", l)
            is_parse_cell = true
        end

        if is_parse_cell && N_parse_cell <= 3
            if N_parse_cell == 0
                N_parse_cell = N_parse_cell + 1
                continue
            end
            ll = split(l, " ", keepempty=false)
            LatVecs[1,N_parse_cell] = parse(Float64, ll[1])
            LatVecs[2,N_parse_cell] = parse(Float64, ll[2])
            LatVecs[3,N_parse_cell] = parse(Float64, ll[3])
            N_parse_cell = N_parse_cell + 1
        end



        #
        # Read atomic positions
        #
        # XXX: Natoms should be read before
        if occursin("ATOMIC_POSITIONS", l)
            is_parse_xyz = true
            if occursin("crystal", l)
                in_fraction = true
                xyz_string_frac = xyz_string_frac*string(Natoms)*"\n\n"
            end
            if occursin("angstrom", l)
                in_angstrom = true
                xyz_string_angstrom = xyz_string_angstrom*string(Natoms)*"\n\n"
            end
            if occursin("bohr", l)
                in_bohr = true
                xyz_string_bohr = xyz_string_bohr*string(Natoms)*"\n\n"
            end
        end

        if is_parse_xyz && (N_parse_xyz <= Natoms) && in_fraction
            if N_parse_xyz == 0
                N_parse_xyz = N_parse_xyz + 1
                continue
            end
            xyz_string_frac = xyz_string_frac*l*"\n"
            N_parse_xyz = N_parse_xyz + 1
        end

        # Atomic positions are given in angstrom
        if is_parse_xyz && (N_parse_xyz <= Natoms) && in_angstrom
            if N_parse_xyz == 0
                N_parse_xyz = N_parse_xyz + 1
                continue
            end
            xyz_string_angstrom = xyz_string_angstrom*l*"\n"
            N_parse_xyz = N_parse_xyz + 1
        end

        # Atomic positions are given in bohr
        if is_parse_xyz && (N_parse_xyz <= Natoms) && in_bohr
            if N_parse_xyz == 0
                N_parse_xyz = N_parse_xyz + 1
                continue
            end
            xyz_string_bohr = xyz_string_bohr*l*"\n"
            N_parse_xyz = N_parse_xyz + 1
        end

        #
        # Read atomic species
        #
        # FIXME: Nspecies must be read before
        if occursin("ATOMIC_SPECIES", l)
            is_parse_species = true
        end

        # Atomic species information (pseudopotentials)
        if is_parse_species && (N_parse_species <= Nspecies)
            if N_parse_species == 0
                N_parse_species = N_parse_species + 1
                continue
            end
            ll = split(l, " ", keepempty=false)
            push!(species_symbols, ll[1])
            push!(species_masses, parse(Float64, ll[2]))
            push!(pspfiles, joinpath(pseudo_dir, ll[3]))
            N_parse_species = N_parse_species + 1
        end

        #
        # Read kpoints
        # Known limitation: only Monkhorst-Pack format is supported
        # or K_POINTS automatic
        #
        if occursin("K_POINTS", l)
            is_parse_kpoints = true
        end

        if is_parse_kpoints && N_parse_kpoints <= 1
            if N_parse_kpoints == 0
                N_parse_kpoints = N_parse_kpoints + 1
                continue
            end
            ll = split(l, " ", keepempty=false)
            meshk1 = parse(Int64, ll[1])
            meshk2 = parse(Int64, ll[2])
            meshk3 = parse(Int64, ll[3])
            N_parse_kpoints = N_parse_kpoints + 1
            #
            # FIXME: Shifts are not read yet
            #
        end

    end
    close(f)


    println(species_symbols)
    println(species_masses)
    println(pspfiles)

    #
    # Some sanity checks
    #

    if ecutwfc <= 0.0
        error("Cannot read ecutwfc, please check or reformat the file")
    end

    if Nspecies <= 0.0
        error("Cannot read Nspecies (ntyp), please check or reformat the file")
    end


    #
    # Set some values if not given using default setting
    #

    # Set ecutrho using default value
    if ecutrho <= 0.0
        ecutrho = 4*ecutwfc
    end

    # Only scale LatVecs with acell if it is a positive value
    if acell > 0.0
        LatVecs = acell*LatVecs
    end

    if in_fraction
        println(xyz_string_frac)
        atoms = init_atoms_xyz_string( xyz_string_frac, in_bohr=true )
        atoms.positions = LatVecs*atoms.positions # convert to bohr
    #
    elseif in_angstrom
        println(xyz_string_angstrom)
        atoms = init_atoms_xyz_string( xyz_string_angstrom )
    #
    elseif in_bohr
        println(xyz_string_bohr)
        atoms = init_atoms_xyz_string( xyz_string_bohr, in_bohr=true )
    else
        error("Cannot read atomic positions")
    end
    # Set unit lattice vectors manually
    atoms.LatVecs = LatVecs

    # Don't forget to convert ecutwfc and ecutrho to Ha
    return PWSCFInput(
        atoms, 0.5*ecutwfc, 0.5*ecutrho, pspfiles,
        (meshk1, meshk2, meshk3)
    )
end