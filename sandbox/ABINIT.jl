function write_abinit( Ham::Hamiltonian;
                       prefix_dir="./",
                       use_smearing=false, kT=0.001 )

    #TODO: Check if prefix is exist or not
    #TODO: If prefix is exist, check if it is directory or not
    #TODO: Create prefix if it is

    pw = Ham.pw
    atoms = Ham.atoms
    electrons = Ham.electrons
    LatVecs = pw.LatVecs
    kpoints = pw.gvecw.kpoints    
    Nspecies = atoms.Nspecies
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    atpos = atoms.positions

    f = open( prefix_dir*"/FILES", "w" )

    println( f, "INPUT" )
    println( f, "LOG1" )
    println( f, "ABINIT_i" )
    println( f, "ABINIT_o" )
    println( f, "ABINIT_" )
    
    # pseudopotential files
    for isp = 1:Nspecies
        println(f, Ham.atoms.SpeciesSymbols[isp]*".psp10")
        write_psp10(Ham.pspots[isp], prefix=prefix_dir)
    end
    close( f )


    # write INPUT file

    f = open( prefix_dir*"/INPUT", "w" )

    println( f, "acell 1.0 1.0 1.0" )
    println( f, "rprim" )
    for i = 1:3
        @printf( f, "%18.10f %18.10f %18.10f\n", LatVecs[1,i], LatVecs[2,i], LatVecs[3,i] )
    end

    @printf(f, "natom %d\n", Natoms)
    znucl = get_Zatoms( atoms )
    @printf(f, "ntypat %d\n", Nspecies)
    @printf(f, "znucl")
    for isp = 1:Nspecies
        @printf(f, " %d", znucl[isp])
    end
    @printf(f, "\n")

    @printf(f, "typat\n")
    for ia = 1:Natoms
        @printf(f, "%d\n", atm2species[ia])
    end

    @printf(f, "xcart\n")
    for ia = 1:Natoms
        @printf(f, "%18.10f %18.10f %18.10f\n", atpos[1,ia], atpos[2,ia], atpos[3,ia])
    end

    @printf(f, "ecut %f\n", pw.ecutwfc)

    if use_smearing
        println(f, "occopt 3")
        @printf(f, "tsmear %f\n", kT)
    else
        println(f, "occopt 1")
    end

    @printf(f, "nband %d\n", electrons.Nstates)

    @printf(f, "ngkpt %d %d %d\n", kpoints.mesh[1], kpoints.mesh[2], kpoints.mesh[3])
    println(f, "chksymbreak 0")
    println(f, "chkprim 0")
    println(f, "nshiftk 1")
    println(f, "shiftk  0.0  0.0  0.0")

    println(f, "nstep 100")
    println(f, "toldfe 1.0d-6")
    println(f, "diemac 12.0")

    if Ham.xcfunc == "VWN"    
        println(f, "ixc -001007")
    elseif Ham.xcfunc == "PBE"
        println(f, "ixc -101130")
    end

    close(f)
end

function read_abinit_etotal( filename::String )

    E_kin      = 0.0
    E_hartree  = 0.0
    E_xc       = 0.0
    E_ewald    = 0.0
    E_pspCore  = 0.0
    E_Ps_loc   = 0.0
    E_Ps_nloc  = 0.0
    mTS        = 0.0
    E_total    = 0.0

    f = open(filename, "r")
    while !eof(f)
        l = readline(f)
        if occursin("Kinetic energy", l)
            E_kin = parse( Float64, split(l, "=")[2] )
        end
        if occursin("Hartree energy", l)
            E_hartree = parse( Float64, split(l, "=")[2] )
        end
        if occursin("XC energy", l)
            E_xc = parse( Float64, split(l, "=")[2] )
        end        
        if occursin("Ewald energy", l)
            E_ewald = parse( Float64, split(l, "=")[2] )
        end
        if occursin("PspCore energy", l)
            E_pspCore = parse( Float64, split(l, "=")[2] )
        end
        if occursin("Loc. psp. energy", l)
            E_Ps_loc = parse( Float64, split(l, "=")[2] )
        end
        if occursin("NL   psp  energy", l)
            E_Ps_nloc = parse( Float64, split(l, "=")[2] )
        end
        if occursin("-kT*entropy", l)
            mTS = parse( Float64, split(l, "=")[2] )
        end
        if occursin("Etotal", l)
            E_total = parse( Float64, split(l, "=")[2] )
        end
    end
    close(f)

    energies = Energies( E_kin, E_Ps_loc, E_Ps_nloc, E_hartree, E_xc, E_ewald, E_pspCore, mTS )
    return energies
end


