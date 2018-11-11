import InteractiveUtils
using Printf
using Random
using BenchmarkTools
using PWDFT

function time_gen_lattice()

    @printf("\n")
    @printf("------------------------------\n")
    @printf("Timing gen_lattice_* functions\n")
    @printf("------------------------------\n")
    @printf("\n")

    @printf("gen_lattice_sc  : ")
    @btime gen_lattice_sc(10.0)

    @printf("gen_lattice_fcc : ")
    @btime gen_lattice_fcc(10.0)
        
    @printf("gen_lattice_bcc : ")
    @btime gen_lattice_bcc(10.0)
end

function time_Atoms()

    @printf("\n")
    @printf("----------------------\n")
    @printf("Timing Atoms functions\n")
    @printf("----------------------\n")
    @printf("\n")

    @printf("Using default constructor (dummy Atoms)          : ")
    @btime Atoms()

    @printf("Using xyz_string with two atoms, default LatVecs : ")
    @btime Atoms(xyz_string="""
    2

    H   0.0   0.0   0.0
    Cl  1.0   0.0   0.0
    """)

    @printf("Using xyz_string with two atoms, given LatVecs   : ")
    @btime Atoms(xyz_string="""
    2

    H   0.0   0.0   0.0
    Cl  1.0   0.0   0.0
    """, LatVecs=gen_lattice_sc(10.0))
end


function time_PWGrid()

    @printf("\n")
    @printf("----------------------\n")
    @printf("Timing PWGrid function\n")
    @printf("----------------------\n")
    @printf("\n")

    ECUTWFC = range(15.0, stop=50.0, step=5.0)    
    LatVecs = gen_lattice_sc(16.0)

    for ecutwfc in ECUTWFC
        @printf("ecutwfc = %4.1f Ha : ", ecutwfc)
        @btime PWGrid($ecutwfc, $LatVecs)
    end
end

function time_Hamiltonian()
    
    @printf("\n")
    @printf("---------------------------\n")
    @printf("Timing Hamiltonian function\n")
    @printf("---------------------------\n")
    @printf("\n")

    atoms = Atoms( xyz_string="""
        1
    
        H   0.0   0.0   0.0
        """,
        LatVecs=gen_lattice_sc(16.0) )
    
    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth"]    
    ECUTWFC = range(15.0, stop=50.0, step=5.0)    

    for ecutwfc in ECUTWFC
        @printf("ecutwfc = %4.1f Ha : ", ecutwfc)
        @btime Hamiltonian( $atoms, $pspfiles, $ecutwfc )
    end

end


function size_Hamiltonian()
    
    @printf("\n")
    @printf("---------------------\n")
    @printf("Memory of Hamiltonian\n")
    @printf("---------------------\n")
    @printf("\n")

    atoms = Atoms( xyz_string="""
        1
    
        H   0.0   0.0   0.0
        """,
        LatVecs=gen_lattice_sc(16.0) )
    
    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth"]    
    ECUTWFC = range(15.0, stop=50.0, step=5.0)    

    for ecutwfc in ECUTWFC
        mem_MB = Base.summarysize( Hamiltonian( atoms, pspfiles, ecutwfc ) )/1024/1024
        @printf("ecutwfc = %4.1f Ha : %f MB\n", ecutwfc, mem_MB)
    end

end

function time_op_K()
    atoms = Atoms(xyz_string_frac=
        """
        1

        Pt  0.0  0.0  0.0
        """, LatVecs=gen_lattice_fcc(3.9231*ANG2BOHR))
    
    pspfiles = ["../pseudopotentials/pade_gth/Pt-q18.gth"]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc,
                       meshk=[8,8,8], extra_states=4 )
        
    # Shortcuts
    pw = Ham.pw
    Focc = Ham.electrons.Focc
    CellVolume = pw.CellVolume
    Ns = pw.Ns

    Random.seed!(4321)
    psiks = rand_BlochWavefunc( Ham )
    
    Rhoe = calc_rhoe( Ham, psiks )
    update!( Ham, Rhoe )

    @printf("Integ rhoe = %18.10f\n\n", sum(Rhoe)*CellVolume/prod(Ns))
    
    @printf("Kinetic operator              : ")
    @btime op_K( $Ham, $psiks )

    @printf("V_Ps_loc                      : ")
    @btime op_V_Ps_loc( $Ham, $psiks )

    @printf("V_loc (Ps loc + Hartree + XC) : ")
    @btime op_V_loc( $Ham, $psiks )

    @printf("V_Ps_nloc                     : ")
    @btime op_V_Ps_nloc( $Ham, $psiks )

    @printf("Hamiltonian                   : ")
    @btime op_H( $Ham, $psiks )


end

println()
InteractiveUtils.versioninfo()

time_gen_lattice()
time_Atoms()
time_PWGrid()
time_Hamiltonian()
size_Hamiltonian()
time_op_K()
