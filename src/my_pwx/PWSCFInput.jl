# using Logging
# It seems that @info and other Logging functions does not play nicely with println?
# The output also is not written to file when using tee

struct PWSCFInput
    atoms::Atoms
    ecutwfc::Float64
    ecutrho::Float64
    pspfiles::Vector{String}
    meshk::Tuple{Int64,Int64,Int64}
    nbnd::Int64
    occupations::String
    smearing::String
    degauss::Float64
    input_dft::String
    nr1::Int64
    nr2::Int64
    nr3::Int64
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

    cell_in_angstrom = false

    is_parse_species = false
    N_parse_species = 0
    pspfiles = String[]
    species_symbols = String[]
    species_masses = Float64[]
    pseudo_dir = "./"

    ibrav = 0

    nbnd = -1 # invalid value

    occupations = ""
    smearing = ""
    degauss = 0.0

    celldm_1 = -1.0
    input_dft = ""

    IGNORED = ["", "\n", "/", "&ELECTRONS", "&CONTROL", "&SYSTEM"]

    # Default values in PWDFT.jl
    # If any of nr1, nr2, nr3 is zero then the FFT grid size
    # will be determined automatically
    nr1 = 0
    nr2 = 0
    nr3 = 0


    f = open(filename, "r")
    
    while !eof(f)
        
        l = readline(f)
        
        # FIXME: This is not robust
        if occursin("  A =", l)
            ll = split(l, "=")
            acell = parse(Float64,ll[end])*ANG2BOHR
            println("Read A = ", acell)
            continue
        end

        # Try reading celldm(1)
        if occursin("  celldm(1) =", l)
            ll = split(l, "=")
            celldm_1 = parse(Float64,ll[end])
            println("Read celldm(1) = ", celldm_1)
            continue
        end

        # TODO: read bcell, ccell, cosAB, cosBC, cosAC
        # TODO: read celldm

        # Read Bravais lattice info (ibrav)
        if occursin(" ibrav =", l)
            ll = split(l, "=", keepempty=false)
            ibrav = parse(Int64, ll[end])
            println("Read ibrav = ", ibrav)
            continue
        end

        # Read number of atoms
        if occursin("nat =", l)
            ll = split(l, "=", keepempty=false)
            Natoms = parse(Int64, ll[end])
            println("Read Natoms = ", Natoms)
            continue
        end

        # Read number of species (ntyp)
        if occursin("ntyp =", l)
            ll = split(l, "=", keepempty=false)
            Nspecies = parse(Int64, ll[end])
            println("Read Nspecies = ", Nspecies)
            continue
        end


        # Read ecutwfc
        if occursin("ecutwfc =", l)
            ll = split(l, "=", keepempty=false)
            ecutwfc = parse(Float64, ll[end])
            println("Read ecutwfc (in Ry) = ", ecutwfc)
            continue
        end

        # Read ecutrho
        if occursin("ecutrho =", l)
            ll = split(l, "=", keepempty=false)
            ecutrho = parse(Float64, ll[end])
            println("Read ecutrho (in Ry) = ", ecutrho)
            continue
        end

        # Occupations
        if occursin("occupations =", l)
            ll = split(l, "=", keepempty=false)
            occupations = replace(replace(ll[end], "'" => ""), "\"" => "") |> strip
            # Be careful of space, need to strip the string
            println("occupations = ", occupations)
            continue
        end

        # degauss
        if occursin("degauss =", l)
            ll = split(l, "=", keepempty=false)
            degauss = parse(Float64, ll[end])
            println("degauss = ", degauss)
            continue
        end

        if occursin("nbnd =", l)
            ll = split(l, "=", keepempty=false)
            nbnd = parse(Int64, ll[end])
            println("nbnd = ", nbnd)
            continue
        end

        # TODO: reduce redundant codes
        if occursin("nr1 =", l)
            ll = split(l, "=", keepempty=false)
            nr1 = parse(Int64, ll[end])
            println("nr1 = ", nbnd)
            continue
        end

        if occursin("nr2 =", l)
            ll = split(l, "=", keepempty=false)
            nr2 = parse(Int64, ll[end])
            println("nr2 = ", nbnd)
            continue
        end

        if occursin("nr3 =", l)
            ll = split(l, "=", keepempty=false)
            nr3 = parse(Int64, ll[end])
            println("nr3 = ", nbnd)
            continue
        end


        # Read pseudo_dir
        if occursin("pseudo_dir =", l)
            # FIXME: Notice the spaces around the equal sign
            # This is to prevent space before / in the path string
            ll = split(l, " = ", keepempty=false)
            pseudo_dir = replace(replace(ll[end], "'" => ""), "\"" => "")
            println("Read pseudo_dir = ", pseudo_dir)
            continue
        end

        # input_dft
        if occursin("input_dft =", l)
            ll = split(l, " = ", keepempty=false)
            input_dft = replace(replace(ll[end], "'" => ""), "\"" => "")
            println("Read input_dft = ", input_dft)
            continue
        end


        #
        # Read cell parameters
        #
        # FIXME: They are assumed to be given in bohr or angstrom!!!!
        if occursin("CELL_PARAMETERS", l)
            is_parse_cell = true
            if occursin("angstrom", l)
                cell_in_angstrom = true
            end
            if occursin("bohr", l)
                cell_in_angstrom = false
            end
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
            continue
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
            continue
        end

        # Atomic positions are given in angstrom
        if is_parse_xyz && (N_parse_xyz <= Natoms) && in_angstrom
            if N_parse_xyz == 0
                N_parse_xyz = N_parse_xyz + 1
                continue
            end
            xyz_string_angstrom = xyz_string_angstrom*l*"\n"
            N_parse_xyz = N_parse_xyz + 1
            continue
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
        # FIXME: The order of pseudopotentials must be the same as in atomic coordinates
        # Otherwise, the  pseudopotentials might get mixed up.
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
            continue
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
            #
            continue
        end

        if !(l in IGNORED)
            println("INFO: Line = $(l) is not processed")
        end

    end
    close(f)

    # Convert LatVecs to bohr if required
    if cell_in_angstrom
        LatVecs[:,:] = LatVecs[:,:]*ANG2BOHR
    end


    println("species_symbols = ", species_symbols)
    println("species_masses = ", species_masses)
    println("pspfiles = ", pspfiles)

    #
    # Some sanity checks
    #

    if ecutwfc <= 0.0
        println("-------------------------------------------------------")
        println("Cannot read ecutwfc, please check or reformat the file")
        println("-------------------------------------------------------")
        error()
    end

    if Nspecies <= 0.0
        println("-------------------------------------------------------")
        println("Cannot read Nspecies (ntyp), please check or reformat the file")
        println("-------------------------------------------------------")
        error()
    end

    if (ibrav == 0) && (!is_parse_cell)
        println("-------------------------------------------------------")
        println("acell (A) must be provided because ibrav = ", ibrav)
        println("-------------------------------------------------------")
        error()
    end

    # TODO: also check for other parameters: esp. celldm, etc
    # Some ibrav need other parameters beside acell
    if (ibrav != 0) && (acell <= 0.0) && (celldm_1 <= 0)
        println("-------------------------------------------------------")
        println("CELL_PARAMETERS must be provided because ibrav = ", ibrav)
        println("-------------------------------------------------------")
        error()
    end 



    #
    # Set some values if not given using default setting
    #

    # Set ecutrho using default value
    if ecutrho <= 0.0
        ecutrho = 4*ecutwfc
    end

    #
    # Prepare for Atoms
    #

    # Only scale LatVecs with acell if it is a positive value
    if ( acell > 0.0 ) && is_parse_cell
        LatVecs = acell*LatVecs
    end

    # FIXME: Special case for commonly used ibrav
    # TODO: Need better code for these 
    if ibrav == 2
        println("INFO: Generating FCC lattice vectors")
        if acell > 0.0
            LatVecs = gen_lattice_fcc(acell)
        elseif celldm_1 > 0.0
            LatVecs = gen_lattice_fcc(celldm_1)
        else
            println("-------------------------------------------------------")
            println("ibrav = ", ibrav, " but both acell and celldm_1 is are not valid")
            println("acell = ", acell, " ibrav = ", ibrav)
            println("-------------------------------------------------------")
            error()
        end
    elseif ibrav == 1
        println("INFO: Generating SC lattice vectors")
        if acell > 0.0
            LatVecs = gen_lattice_sc(acell)
        elseif celldm_1 > 0.0
            LatVecs = gen_lattice_sc(celldm_1)
        else
            println("-------------------------------------------------------")
            println("ibrav = ", ibrav, " but both acell and celldm_1 is are not valid")
            println("acell = ", acell, " ibrav = ", ibrav)
            println("-------------------------------------------------------")
            error()
        end
    else
        @info("\nPlease check if the following LatVecs makes sense:\n")
        @info(LatVecs)
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

    println(atoms)

    # Don't forget to convert ecutwfc and ecutrho to Ha
    return PWSCFInput(
        atoms, 0.5*ecutwfc, 0.5*ecutrho, pspfiles,
        (meshk1, meshk2, meshk3),
        nbnd, occupations, smearing, degauss, input_dft,
        nr1, nr2, nr3
    )
end