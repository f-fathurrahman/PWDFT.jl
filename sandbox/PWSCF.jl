"""
Write GTH pseudopotential in the format that can be read by PWSCF.
"""
function write_psp_pwscf( psp::PsPot_GTH; prefix="." )
    atsymb = psp.atsymb
    zval = psp.zval
    zatom = PWDFT.ZATOMS[atsymb]
    lmax = psp.lmax
    
    f = open(prefix*"/"*atsymb*".gth", "w")
    @printf(f, "Hartwigsen-Goedecker-Hutter psp for %s, from PRB58, 3641 (1998)\n", atsymb)
    @printf(f, "%d %d 010605\n", Int(zatom), Int(zval))
    @printf(f, "10 1 %d 0 2001 0\n", lmax)
    
    # local pseudopotential
    rlocal = psp.rlocal
    c = psp.c
    ncoef = count(abs.(c) .> 0.0)
    @printf(f, "%18.10f %5d", rlocal, ncoef)
    for i = 1:ncoef
        @printf(f, "%18.10f ", c[i])
    end
    @printf(f, "\n")


    # Nonlocal part
    Nproj_l = psp.Nproj_l
    rc = psp.rc
    h = psp.h
    @printf(f, "%5d\n", lmax+1)
    for l = 0:lmax
        @printf(f, "%18.10f %5d", rc[l+1], Nproj_l[l+1])
        for i = 1:Nproj_l[l+1]
            for j = i:Nproj_l[l+1]
                @printf(f, "%18.10f ", h[l+1,i,j])
            end
            @printf(f, "               \n")
        end
        if l > 0
            for i = 1:Nproj_l[l+1]
                for j = i:Nproj_l[l+1]
                    @printf(f, "%18.10f ", 0.0)
                end
                @printf(f, "               \n")
            end
        end
    end

    close(f)

    return
end

"""
Write a PWSCF input file from a given `Hamiltonian` instance.
Limitation:
- Calculation is `scf`
- Atomic mass is set to 1.0
"""
function write_pwscf( Ham::Hamiltonian; filename="PWINPUT",
                      prefix="./", prefix_psp="./",
                      use_smearing=false, kT=0.001 )
    atoms = Ham.atoms
    pw = Ham.pw

    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    SpeciesSymbols = atoms.SpeciesSymbols
    atpos = atoms.positions
    atsymbs = atoms.atsymbs

    ecutwfc = pw.ecutwfc

    f = open(prefix*filename, "w")

    @printf(f, "&CONTROL\n")
    @printf(f, "  calculation = 'scf'\n")
    @printf(f, "  restart_mode = 'from_scratch'\n")
    @printf(f, "  pseudo_dir = './'\n")
    @printf(f, "  outdir = './tmp'\n")
    @printf(f, "  verbosity = 'high'\n")
    @printf(f, "  disk_io = 'none'\n")
    @printf(f, "/\n\n")

    @printf(f, "&SYSTEM\n")
    @printf(f, "  ibrav = 0\n")
    @printf(f, "  nat = %d\n", Natoms)
    @printf(f, "  ntyp = %d\n", Nspecies)
    @printf(f, "  ecutwfc = %18.10f\n", ecutwfc*2)
    if use_smearing
        @printf(f, "  occupations = 'smearing'\n")
        @printf(f, "  smearing = 'fermi-dirac'\n")
        @printf(f, "  degauss = %18.10f\n", kT*2)
    end
    @printf(f, "/\n\n")

    @printf(f, "&ELECTRONS\n")
    @printf(f, "  electron_maxstep = 150\n")
    @printf(f, "  mixing_beta = 0.1\n")
    @printf(f, "/\n\n")

    @printf(f, "ATOMIC_SPECIES\n")
    for isp = 1:Nspecies
        ss = SpeciesSymbols[isp]
        @printf(f, "%5s 1.0 %s\n", ss, ss*".gth")
    end
    @printf(f, "\n")

    @printf(f, "ATOMIC_POSITIONS bohr\n")
    for ia = 1:Natoms
        @printf(f, "%s %18.10f %18.10f %18.10f\n",
                atsymbs[ia], atpos[1,ia], atpos[2,ia], atpos[3,ia])
    end
    @printf(f, "\n")

    kmesh = pw.gvecw.kpoints.mesh
    @printf(f, "K_POINTS automatic\n")
    @printf(f, "%d %d %d 0 0 0\n", kmesh[1], kmesh[2], kmesh[3])
    @printf(f, "\n")

    LatVecs = pw.LatVecs
    @printf(f, "CELL_PARAMETERS bohr\n")
    for i = 1:3
        @printf(f, "%18.10f %18.10f %18.10f\n", LatVecs[1,i], LatVecs[2,i], LatVecs[3,i])
    end
    @printf(f, "\n")

    close(f)

    for psp in Ham.pspots
        write_psp_pwscf(psp, prefix=prefix_psp)
    end

end

mutable struct EnergiesPWSCF
    OneElectron::Float64
    Hartree::Float64
    XC::Float64
    NN::Float64
    mTS::Float64
end

import Base: println
function println( energies::EnergiesPWSCF; use_smearing=false )

    @printf("one-elec   energy: %18.10f\n", energies.OneElectron )
    @printf("Hartree    energy: %18.10f\n", energies.Hartree )
    @printf("XC         energy: %18.10f\n", energies.XC )

    if ( abs(energies.mTS) > eps() ) || use_smearing
        @printf("-TS              : %18.10f\n", energies.mTS)
    end

    @printf("-------------------------------------\n")
    
    E_elec = energies.OneElectron + energies.Hartree + energies.XC + energies.mTS
    
    @printf("Electronic energy: %18.10f\n", E_elec)
    @printf("NN         energy: %18.10f\n", energies.NN )
    @printf("-------------------------------------\n")
    
    E_total = E_elec + energies.NN
    
    if use_smearing
        @printf("Total free energy: %18.10f\n", E_total)
    else
        @printf("Total      energy: %18.10f\n", E_total )
    end
end    



function read_pwscf_etotal( filename::String )
    E_one_ele = 0.0
    E_hartree = 0.0
    E_xc = 0.0
    E_nn = 0.0
    mTS = 0.0

    # XXX: Don't forget to convert to Hartree

    f = open(filename, "r")
    while !eof(f)
        l = readline(f)
        if occursin("one-electron contribution", l)
            ll = split(l, " ", keepempty=false)
            E_one_ele = parse( Float64, ll[4] )*0.5
        end
        if occursin("hartree contribution", l)
            ll = split(l, " ", keepempty=false)
            E_hartree = parse( Float64, ll[4] )*0.5
        end
        if occursin("xc contribution", l)
            ll = split(l, " ", keepempty=false)
            E_xc = parse( Float64, ll[4] )*0.5
        end
        if occursin("ewald contribution", l)
            ll = split(l, " ", keepempty=false)
            E_nn = parse( Float64, ll[4] )*0.5
        end        
        if occursin("smearing contrib", l)
            ll = split(l, " ", keepempty=false)
            mTS = parse( Float64, ll[5] )*0.5
        end
    end

    return EnergiesPWSCF(E_one_ele, E_hartree, E_xc, E_nn, mTS)
end




"""
Parse PWSCF input file generated by cif2cell program.
Returns an instance of `Atoms` and KPoints' mesh.
"""
function read_pwscf_cif2cell( filename::String )

    # Default values, some are not valid for PWSCF
    acell = -1.0

    Natoms = 0

    LatVecs = zeros(3,3)
    is_parse_cell = false
    N_parse_cell = 0

    xyz_string_frac = ""
    is_parse_xyz = false
    N_parse_xyz = 0

    is_parse_kpoints = false
    N_parse_kpoints = 0
    meshk = [0, 0, 0]

    f = open(filename, "r")
    
    while !eof(f)
        
        l = readline(f)
        println()
        println("Current line = ", l)
        
        if occursin("  A =", l)
            ll = split(l, "=")
            println(ll)
            acell = parse(Float64,ll[end])*ANG2BOHR
            println("acell = ", acell)
        end

        if occursin("nat =", l)
            ll = split(l, "=", keepempty=false)
            Natoms = parse(Int64, ll[end])
            xyz_string_frac = xyz_string_frac*string(Natoms)*"\n\n"
        end

        if occursin("CELL_PARAMETERS", l)
            is_parse_cell = true
        end

        if is_parse_cell && N_parse_cell <= 3
            if N_parse_cell == 0
                N_parse_cell = N_parse_cell + 1
                continue
            end
            println("Parsing cell")
            ll = split(l, " ", keepempty=false)
            println(ll)
            LatVecs[1,N_parse_cell] = parse(Float64, ll[1])
            LatVecs[2,N_parse_cell] = parse(Float64, ll[2])
            LatVecs[3,N_parse_cell] = parse(Float64, ll[3])
            N_parse_cell = N_parse_cell + 1
        end


        if occursin("ATOMIC_POSITIONS", l)
            is_parse_xyz = true
        end

        if is_parse_xyz && N_parse_xyz <= Natoms
            if N_parse_xyz == 0
                N_parse_xyz = N_parse_xyz + 1
                continue
            end
            println("Parsing xyz")
            xyz_string_frac = xyz_string_frac*l*"\n"
            N_parse_xyz = N_parse_xyz + 1
        end

        if occursin("K_POINTS", l)
            is_parse_kpoints = true
        end

        if is_parse_kpoints && N_parse_kpoints <= 1
            if N_parse_kpoints == 0
                N_parse_kpoints = N_parse_kpoints + 1
                continue
            end
            println("Parsing kpoints")
            ll = split(l, " ", keepempty=false)
            meshk[1] = parse(Int64,ll[1])
            meshk[2] = parse(Int64,ll[2])
            meshk[3] = parse(Int64,ll[3])
            println("meshk = ", meshk)
            N_parse_kpoints = N_parse_kpoints + 1
        end

    end
    close(f)

    LatVecs = acell*LatVecs
    print_matrix(LatVecs)
    println(xyz_string_frac)

    atoms = init_atoms_xyz_string( xyz_string_frac, in_bohr=true )
    atoms.positions = LatVecs*atoms.positions
    atoms.LatVecs = LatVecs

    return atoms, meshk
end

