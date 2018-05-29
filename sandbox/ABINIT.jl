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
    psp_path = Array{String,1}(Nspecies)
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

function write_abinit( Ham::PWHamiltonian;
                       abinit_psp=nothing, prefix="../TEMP/", prefix_psp="./" )

    f = open( prefix*"FILES", "w" )

    println( f, "INPUT" )
    println( f, "OUTLOG" )
    println( f, "ABINIT_i" )
    println( f, "ABINIT_o" )
    println( f, "ABINIT_" )

    Nspecies = Ham.atoms.Nspecies
    if abinit_psp == nothing
        abinit_psp = get_abinit_psp(Ham.atoms.SpeciesSymbols)
    end
    println(abinit_psp)
    @assert length(abinit_psp) == Nspecies
    for isp = 1:Nspecies
        println( f, prefix_psp*abinit_psp[isp] )
    end
    close( f )

    f = open( prefix*"INPUT", "w" )
    println( f, "acell 1.0 1.0 1.0" )
    println( f, "rprim" )
    LatVecs = Ham.pw.LatVecs
    for i = 1:3
        @printf( f, "%18.10f %18.10f %18.10f\n", LatVecs[1,i], LatVecs[2,i], LatVecs[3,i] )
    end
    close(f)
end


function test_main()
    # initialize atoms and Hamiltonian
    atoms = init_atoms_xyz("../structures/CuSO4.xyz")
    ecutwfc = 15.0
    Ham = PWHamiltonian( atoms, ecutwfc )

    write_abinit( Ham )
end

test_main()
