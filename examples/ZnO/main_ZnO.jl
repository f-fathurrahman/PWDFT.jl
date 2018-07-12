using Printf
using PWDFT

function test_main( ; method="SCF" )
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        4

        Zn      0.3333333   0.6666667   0.0000000
        Zn      0.6666667   0.3333333   0.5000000
        O       0.3333333   0.6666667   0.3450000
        O       0.6666667   0.3333333   0.8450000
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_hexagonal( 3.2495*ANG2BOHR, 5.2069*ANG2BOHR )
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)
    write_xsf( "TEMP_ZnO.xsf", atoms )

    # Initialize Hamiltonian
    pspfiles = ["../../pseudopotentials/pade_gth/Zn-q2.gth",
                "../../pseudopotentials/pade_gth/O-q6.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, meshk=[3,3,3], verbose=true )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

    #
    # Solve the KS problem
    #
    if method == "SCF"
        KS_solve_SCF!( Ham )

    elseif method == "Emin"
        KS_solve_Emin_PCG!( Ham, verbose=true )

    elseif method == "DCM"
        KS_solve_DCM!( Ham, NiterMax=15 )

    else
        println("ERROR: unknow method = ", method)
    end

    println("\nTotal energy components")
    println(Ham.energies)

end

@time test_main(method="Emin")
@time test_main(method="SCF")
@time test_main(method="DCM")


#=
!    total energy              =    -234.44638420 Ry = -117.2231921 Ha
     one-electron contribution =     -22.62900357 Ry
     hartree contribution      =      43.45832362 Ry
     xc contribution           =     -38.16290548 Ry
     ewald contribution        =    -217.11279878 Ry = -108.55639939 Ha


ABINIT Zn-q12
    Kinetic energy  =  7.30418778115223E+01
    Hartree energy  =  2.17431897202879E+01
    XC energy       = -1.90781737938864E+01
    Ewald energy    = -1.08556401221583E+02
    PspCore energy  =  4.40899664450704E+00
    Loc. psp. energy= -8.81088473603044E+01
    NL   psp  energy= -6.67842381110405E-01
    >>>>>>>>> Etotal= -1.17217200580567E+02

ABINIT Zn-q2
    Kinetic energy  =  1.92461355181398E+01
    Hartree energy  =  9.21163371686482E+00
    XC energy       = -7.45191903004917E+00
    Ewald energy    = -2.34051747159905E+01
    PspCore energy  =  4.13113070747948E-01
    Loc. psp. energy= -3.40424483141658E+01
    NL   psp  energy=  3.53153690538100E+00
    >>>>>>>>> Etotal= -3.24971228490719E+01

With Zn-q2
Kinetic    energy:      19.2479527554
Ps_loc     energy:     -33.6337949101
Ps_nloc    energy:       3.5303208825
Hartree    energy:       9.2217086060
XC         energy:      -7.4648139530
-------------------------------------
Electronic energy:      -9.0986266191
NN         energy:     -23.4051760908
-------------------------------------
Total      energy:     -32.5038027099


With Zn-q12:
Kinetic    energy:      73.0254762512
Ps_loc     energy:     -83.7155981990
Ps_nloc    energy:      -0.6318222977
Hartree    energy:      21.7499768105
XC         energy:     -19.1016050852
-------------------------------------
Electronic energy:      -8.6735725202
NN         energy:    -108.5564075980
-------------------------------------
Total      energy:    -117.2299801183
=#
