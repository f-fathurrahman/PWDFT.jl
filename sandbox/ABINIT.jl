using Printf
using PWDFT

const ABINIT_PSP_LDA = 
"""
Ag.11.LDA.hgh
Al.3.LDA.hgh
Ar.8.LDA.hgh
As.5.LDA.hgh
At.7.LDA.hgh
Au.11.LDA.hgh
B.3.LDA.hgh
Ba.10.LDA.hgh
Ba.2.LDA.hgh
Be.4.LDA.hgh
Bi.5.LDA.hgh
Br.7.LDA.hgh
C.4.LDA.hgh
Ca.10.LDA.hgh
Cd.12.LDA.hgh
Ce.12.LDA.hgh
Cl.7.LDA.hgh
Co.17.LDA.hgh
Cr.14.LDA.hgh
Cs.9.LDA.hgh
Cu.11.LDA.hgh
Dy.20.LDA.hgh
Er.22.LDA.hgh
Eu.17.LDA.hgh
F.7.LDA.hgh
Fe.16.LDA.hgh
Ga.13.LDA.hgh
Gd.18.LDA.hgh
Ge.4.LDA.hgh
H.1.LDA.hgh
He.2.LDA.hgh
Hf.12.LDA.hgh
Hg.12.LDA.hgh
Ho.21.LDA.hgh
I.7.LDA.hgh
In.13.LDA.hgh
Ir.17.LDA.hgh
K.9.LDA.hgh
Kr.8.LDA.hgh
La.11.LDA.hgh
Li.3.LDA.hgh
Lu.25.LDA.hgh
Mg.10.LDA.hgh
Mn.15.LDA.hgh
Mo.14.LDA.hgh
N.5.LDA.hgh
Na.9.LDA.hgh
Nb.13.LDA.hgh
Nd.14.LDA.hgh
Ne.8.LDA.hgh
Ni.18.LDA.hgh
O.6.LDA.hgh
Os.16.LDA.hgh
P.5.LDA.hgh
Pb.4.LDA.hgh
Pd.18.LDA.hgh
Pm.15.LDA.hgh
Po.6.LDA.hgh
Pr.13.LDA.hgh
Pt.18.LDA.hgh
Rb.9.LDA.hgh
Re.15.LDA.hgh
Rh.17.LDA.hgh
Rn.8.LDA.hgh
Ru.16.LDA.hgh
S.6.LDA.hgh
Sb.5.LDA.hgh
Sc.11.LDA.hgh
Se.6.LDA.hgh
Si.4.LDA.hgh
Sm.16.LDA.hgh
Sn.4.LDA.hgh
Sr.10.LDA.hgh
Ta.13.LDA.hgh
Tb.19.LDA.hgh
Tc.15.LDA.hgh
Te.6.LDA.hgh
Ti.12.LDA.hgh
Tl.13.LDA.hgh
Tm.23.LDA.hgh
V.13.LDA.hgh
W.14.LDA.hgh
Xe.8.LDA.hgh
Y.11.LDA.hgh
Yb.24.LDA.hgh
Zn.12.LDA.hgh
Zr.12.LDA.hgh
"""

function get_psp_atsymb( path_str )
    return split( path_str, "." )[1]
end

function get_psp_chg( path_str )
    return split( path_str, "." )[2]
end

function get_abinit_psp( SpeciesSymbols )
    psps = split( ABINIT_PSP_LDA )
    Nspecies = size(SpeciesSymbols)[1]
    psp_path = Array{String,1}(undef,Nspecies)
    for isp = 1:Nspecies
        s = SpeciesSymbols[isp]
        for ps in psps
            psp_atsymb = get_psp_atsymb( ps )
            if s == psp_atsymb
                println("Found")
                psp_path[isp] = ps
            end
        end
    end
    return psp_path
end

function write_abinit( Ham::Hamiltonian;
                       abinit_psp=nothing, prefix="./", prefix_psp="./",
                       use_smearing=false, kT=0.001 )

    pw = Ham.pw
    atoms = Ham.atoms
    electrons = Ham.electrons
    LatVecs = pw.LatVecs
    kpoints = pw.gvecw.kpoints    
    Nspecies = atoms.Nspecies
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    atpos = atoms.positions

    f = open( prefix*"FILES", "w" )

    println( f, "INPUT" )
    println( f, "LOG1" )
    println( f, "ABINIT_i" )
    println( f, "ABINIT_o" )
    println( f, "ABINIT_" )

    if abinit_psp == nothing
        abinit_psp = get_abinit_psp(atoms.SpeciesSymbols)
    end
    println(abinit_psp)
    @assert length(abinit_psp) == Nspecies
    for isp = 1:Nspecies
        println( f, prefix_psp*abinit_psp[isp] )
    end
    close( f )


    # write INPUT file

    f = open( prefix*"INPUT", "w" )

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
    println(energies)
    return energies
end


function test_CuSO4()
    # initialize atoms and Hamiltonian
    atoms = init_atoms_xyz("../structures/CuSO4.xyz")
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, ecutwfc )

    write_abinit( Ham, prefix_psp="../_compare/abinit_psp/",
                       prefix="../TEMP_ABINIT/" )
end
#test_CuSO4()


function test_H2()
    # initialize atoms and Hamiltonian
    atoms = Atoms( xyz_file="../structures/H2.xyz",
                   LatVecs=gen_lattice_sc(16.0) )
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, ecutwfc )

    write_abinit( Ham, prefix_psp="../_compare/abinit_psp/",
                       prefix="../TEMP_ABINIT/" )
end
#test_H2()

function test_Al_fcc()

    atoms = Atoms( xyz_string_frac=
    """
    1

    Al  0.0  0.0  0.0
    """, in_bohr=true,
    LatVecs = gen_lattice_fcc(7.6525970200) )
    pspfiles = ["../pseudopotentials/pade_gth/Al-q3.gth"]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, xcfunc="LDA",
                       Nspin=1, meshk=[8,8,8], extra_states=4 )

    write_abinit( Ham, use_smearing=true,
                  prefix_psp="../_compare/abinit_psp/",
                  prefix="../TEMP_ABINIT/" )
end
#test_Al_fcc()

function test_read_abinit_etotal()
    read_abinit_etotal("../TEMP_ABINIT/LOG1")
end
#test_read_abinit_etotal()

function test_Pt_fcc()
    # Atoms
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_fcc(3.9231*ANG2BOHR))

    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/Pt-q18.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, xcfunc="PBE",
                        meshk=[3,3,3], extra_states=4 )
    
    write_abinit( Ham, use_smearing=true,
                  prefix_psp="../_compare/abinit_psp/",
                  prefix="../TEMP_ABINIT/" )
end
test_Pt_fcc()