"""
Write a GPAW input script (Python code) from a given `Hamiltonian` instance.

Internal HGH pseudopotentials in GPAW are used.

FIXME:
The keyword `etot_conv_thr` is in the function argument but it is not used.
There are several convergence criteria for GPAW. The energy criteria is defined
per number of electrons.
For the moment I decided to use the default criteria.
"""
function write_gpaw( Ham::Hamiltonian; filename="main.py",
                     prefix_dir="./",
                     use_smearing=false, kT=0.001,
                     etot_conv_thr=1e-6 )

    atoms = Ham.atoms

    if !Base.Filesystem.isdir(prefix_dir)
        Base.Filesystem.mkdir(prefix_dir)
    end

    f = open( joinpath(prefix_dir, filename), "w")

    @printf(f, "from ase import Atoms\n")
    @printf(f, "from gpaw import GPAW, PW, FermiDirac\n")
    @printf(f, "from ase.units import Bohr, Hartree\n")

    # Write out the system name. GPAW should parse them to determine the atomic
    # species present in the system.
    atoms_name = ""
    for ia = 1:atoms.Natoms
        atoms_name = atoms_name * atoms.SpeciesSymbols[atoms.atm2species[ia]]
    end

    # Atomic positions. The corresponding atomic species should be determined
    # from atoms_name.
    atpos = atoms.positions/ANG2BOHR
    @printf(f, "atoms = Atoms('%s', positions=[\n", atoms_name)
    for ia in 1:atoms.Natoms
        @printf(f, "    [%18.10f, %18.10f, %18.10f],\n", atpos[1,ia], atpos[2,ia], atpos[3,ia])
    end
    @printf(f, "    ],\n")

    # This is our defalt
    @printf(f, "    pbc = [True, True, True],\n")

    # Cell vectors
    v = atoms.LatVecs' / ANG2BOHR
    @printf(f, "    cell=[\n")
    for i in 1:3
        @printf(f, "    [%18.10f, %18.10f, %18.10f],\n", v[i,1], v[i,2], v[i,3])
    end
    @printf(f, "    ]\n")
    @printf(f, ")\n")

    # For visualization
    @printf(f, "atoms.write('ATOMS.xsf')\n")

    @printf(f, "ecutwfc = %f*Hartree\n", Ham.pw.ecutwfc)

    # Start building GPAW calculator
    @printf(f, "calc = GPAW( mode=PW(ecutwfc),\n")

    # HGH pseudopotentials.
    @printf(f, "             setups='hgh',\n")
    
    # Default XC functional (via LibXC)
    @printf(f, "             xc='LDA_X+LDA_C_VWN',\n")
    
    # Spin polarized of not
    if Ham.electrons.Nspin == 2
        @printf(f, "             spinpol=True,\n")
    else
        @printf(f, "             spinpol=False,\n")
    end

    # Number of bands
    @printf(f, "             nbands=%d,\n", Ham.electrons.Nstates)

    # Smearing of not
    if use_smearing
        @printf(f, "             occupations=FermiDirac(%f*Hartree),\n", kT)
    else
        @printf(f, "             occupations=FermiDirac(0),\n")
    end

    # 'gamma': True is set for unshifted kpoint mesh
    meshk = Ham.pw.gvecw.kpoints.mesh
    @printf(f, "             kpts={'size': (%d, %d, %d), 'gamma':True},\n", meshk[1], meshk[2], meshk[3])

    # Print out to stdout
    @printf(f, "             txt='-')\n")

    # Set the calculator and calculate the total energy.
    @printf(f, "atoms.set_calculator(calc)\n")
    @printf(f, "e1 = atoms.get_potential_energy()\n")
    # Convert to Hartree
    @printf(f, "print('\\n!!!!!! Total energy      : %%18.10f eV = %%18.10f Hartree\\n' %% (e1,e1/Hartree))\n")

    close(f)
end


# FIXME: need to know the exact meaning of these terms and their correspondence with PWDFT.jl

struct EnergiesGPAW
    Kinetic::Float64
    Potential::Float64
    External::Float64
    XC::Float64
    Entropy::Float64
    Local::Float64
    Free_energy::Float64
    Extrapolated::Float64
end

# No `sum` method for this type


# Smearing information is always printed out, no `use_smearing` keyword.

import Base: show
function show( io::IO, energies::EnergiesGPAW )
    @printf(io, "Kinetic   energy: %18.10f\n", energies.Kinetic)
    @printf(io, "Potential energy: %18.10f\n", energies.Potential)
    @printf(io, "External  energy: %18.10f\n", energies.External)
    @printf(io, "XC        energy: %18.10f\n", energies.XC)
    @printf(io, "Entropy   energy: %18.10f\n", energies.Entropy)
    @printf(io, "Local     energy: %18.10f\n", energies.Local)
    @printf(io, "------------------------------------\n")
    @printf(io, "Free energy     : %18.10f\n", energies.Free_energy)
    @printf(io, "Extrapolated    : %18.10f\n", energies.Extrapolated)
end

show( energies::EnergiesGPAW ) = show( stdout, energies )

import Base: println
println( energies::EnergiesGPAW ) = show( energies )



function read_gpaw_etotal( filename::String )
    
    Kinetic = 0.0
    Potential = 0.0
    External = 0.0
    XC = 0.0
    Entropy = 0.0 # mTS (ie -TS)
    Local = 0.0
    Free_energy = 0.0
    Extrapolated = 0.0

    f = open(filename, "r")
    while !eof(f)
        
        l = readline(f)
        #
        if occursin("Energy contributions", l)
            #
            l = readline(f) # skip
            #            
            l = readline(f)
            Kinetic = parse( Float64, split(l, keepempty=false)[end] )/(2*Ry2eV)
            #
            l = readline(f)
            Potential = parse( Float64, split(l, keepempty=false)[end] )/(2*Ry2eV)
            #
            l = readline(f)
            External = parse( Float64, split(l, keepempty=false)[end] )/(2*Ry2eV)
            #
            l = readline(f)
            XC = parse( Float64, split(l, keepempty=false)[end] )/(2*Ry2eV)
            #
            l = readline(f)
            Entropy = parse( Float64, split(l, keepempty=false)[end] )/(2*Ry2eV)
            #
            l = readline(f)
            Local = parse( Float64, split(l, keepempty=false)[end] )/(2*Ry2eV)
            #
            l = readline(f)
            #
            l = readline(f)
            Free_energy = parse( Float64, split(l, keepempty=false)[end] )/(2*Ry2eV)
            #
            l = readline(f)
            Extrapolated = parse( Float64, split(l, keepempty=false)[end] )/(2*Ry2eV)
            #
            break
        end

    end

    close(f)
    
    return EnergiesGPAW( Kinetic, Potential, External, XC, Entropy, Local, Free_energy, Extrapolated )
end
