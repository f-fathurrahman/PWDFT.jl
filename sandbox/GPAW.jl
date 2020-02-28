"""
Write a GPAW input script (Python code) from a given `Hamiltonian` instance.

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

    atoms_name = ""
    for ia = 1:atoms.Natoms
        atoms_name = atoms_name * atoms.SpeciesSymbols[atoms.atm2species[ia]]
    end

    atpos = atoms.positions/ANG2BOHR
    @printf(f, "atoms = Atoms('%s', positions=[\n", atoms_name)
    for ia in 1:atoms.Natoms
        @printf(f, "    [%18.10f, %18.10f, %18.10f],\n", atpos[1,ia], atpos[2,ia], atpos[3,ia])
    end
    @printf(f, "    ],\n")

    @printf(f, "    pbc = [True, True, True],\n")

    v = atoms.LatVecs' / ANG2BOHR
    @printf(f, "    cell=[\n")
    for i in 1:3
        @printf(f, "    [%18.10f, %18.10f, %18.10f],\n", v[i,1], v[i,2], v[i,3])
    end
    @printf(f, "    ]\n")
    @printf(f, ")\n")

    @printf(f, "atoms.write('ATOMS.xsf')\n") # for comparison

    @printf(f, "ecutwfc = %f*Hartree\n", Ham.pw.ecutwfc)
    @printf(f, "calc = GPAW( mode=PW(ecutwfc),\n")
    @printf(f, "             setups='hgh',\n")
    @printf(f, "             xc='LDA_X+LDA_C_VWN',\n")
    if Ham.electrons.Nspin == 2
        @printf(f, "             spinpol=True,\n")
    else
        @printf(f, "             spinpol=False,\n")
    end
    if use_smearing
        @printf(f, "             occupations=FermiDirac(%f*Hartree),\n", kT)
    else
        @printf(f, "             occupations=FermiDirac(0),\n")
    end

    # 'gamma': True is set for unshifted kpoint mesh
    meshk = Ham.pw.gvecw.kpoints.mesh
    @printf(f, "             kpts={'size': (%d, %d, %d), 'gamma':True},\n", meshk[1], meshk[2], meshk[3])
    @printf(f, "             txt='-')\n")
    @printf(f, "atoms.set_calculator(calc)\n")
    @printf(f, "e1 = atoms.get_potential_energy()\n")
    @printf(f, "print('\\n!!!!!! Total energy      : %%18.10f eV = %%18.10f Hartree\\n' %% (e1,e1/Hartree))\n")

    close(f)
end