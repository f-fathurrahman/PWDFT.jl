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
function write_pwscf( Ham::Hamiltonian;
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

    f = open(prefix*"PWINPUT", "w")

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

