using Printf
using LinearAlgebra
using Random
using PWDFT

include("dump_bandstructure.jl")
include("gen_kpath.jl")

function test_graphene()

    # Atoms
    atoms = init_atoms_xyz_string(
        """
        2

        C  0.047580  0.024180  0.0
        C  0.047662  1.446213  0.0
        """)
    LatVecs = zeros(3,3)
    LatVecs[:,1] = [2.4638581459, 0.0000000000, 0.0000000000]
    LatVecs[:,2] = [-1.2319290729, 2.1337637457, 0.0000000000]
    LatVecs[:,3] = [0.0, 0.0, 5.0]
    atoms.LatVecs = LatVecs*ANG2BOHR
    println(atoms)
    write_xsf("TEMP_graphene.xsf", atoms)

    pspfiles = ["../pseudopotentials/pade_gth/C-q4.gth"]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc,
                       meshk=[4,4,1], extra_states=1 )
    println(Ham)

    KS_solve_SCF!( Ham )

    #
    # Band structure calculation
    #
    #kpoints, kpt_spec, kpt_spec_labels = kpath_from_file(atoms, "KPATH_GRAPHENE_40")
    kpoints, kpt_spec, kpt_spec_labels = gen_kpath(atoms, "GMKG", "hexagonal")

    # New pw, ecutwfc should be the same as SCF calc
    Ham.pw = PWGrid(ecutwfc, atoms.LatVecs, kpoints=kpoints)

    Ngw = Ham.pw.gvecw.Ngw
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    k = Ham.pw.gvecw.kpoints.k

    # Manually construct Ham.electrons
    Ham.electrons = Electrons( atoms, Ham.pspots, Nspin=1, Nkpt=kpoints.Nkpt,
                               Nstates_empty=1 )

    Nstates = Ham.electrons.Nstates
    ebands = Ham.electrons.ebands
    Nspin = Ham.electrons.Nspin

    Nkspin = Nkpt*Nspin

    # Also, new PsPotNL
    Ham.pspotNL = PsPotNL( Ham.pw, atoms, Ham.pspots, kpoints, check_norm=false )

    psiks = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    evals = zeros(Float64,Nstates,Nkspin)

    Random.seed!(1234)
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        ikspin = ik + (ispin - 1)*Nkpt
        psiks[ikspin] = ortho_sqrt(rand(ComplexF64,Ngw[ik],Nstates))
    end
    end

    k = Ham.pw.gvecw.kpoints.k
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        #
        @printf("\nispin = %d, ik = %d, ikspin=%d, Ngw = %d\n", ispin, ik, ikspin, Ngw[ik])
        @printf("kpts = [%f,%f,%f]\n", k[1,ik], k[2,ik], k[3,ik])
        evals[:,ikspin], psiks[ikspin] =
        diag_LOBPCG( Ham, psiks[ikspin], verbose_last=true )
    end
    end

    dump_bandstructure( evals, kpoints.k, kpt_spec, kpt_spec_labels, filename="TEMP_GRAPHENE.dat" )

end

test_graphene()

