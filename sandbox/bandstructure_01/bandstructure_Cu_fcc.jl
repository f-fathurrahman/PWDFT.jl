using Printf
using LinearAlgebra
using Random
using PWDFT

include("dump_bandstructure.jl")
include("gen_kpath.jl")

function test_Cu_fcc()

    # Atoms
    atoms = init_atoms_xyz_string(
        """
        1

        Cu  0.0  0.0  0.0
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(3.61496*ANG2BOHR)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)

    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/Cu-q11.gth"]
    ecutwfc = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc,
                       meshk=[16,16,16], extra_states=1 )
    println(Ham)

    #
    # Solve the KS problem
    #
    KS_solve_SCF!( Ham, betamix=0.2, mix_method="rpulay", use_smearing=true, kT=0.01 )


    #
    # Band structure calculation
    #
    #kpoints, kpt_spec, kpt_spec_labels = kpath_from_file(atoms, "KPATH_FCC_60")
    kpoints, kpt_spec, kpt_spec_labels = gen_kpath(atoms, "WLGXWK", "fcc", Δk=0.05 )
    #kpoints, kpt_spec, kpt_spec_labels = gen_kpath(atoms, "GXWKGLUWLK", "fcc", Δk=0.05 )

    # New pw
    pw = PWGrid(ecutwfc, atoms.LatVecs, kpoints=kpoints)
    Ham.pw = pw

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
    Ham.pspotNL = PsPotNL( pw, atoms, Ham.pspots, kpoints, check_norm=false )


    psiks = BlochWavefunc(undef,Nkspin)
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
        #evals[:,ikspin], psiks[ikspin] =
        #diag_davidson( Ham, psiks[ikspin], verbose_last=true )
    end
    end

    dump_bandstructure( evals, kpoints.k, kpt_spec, kpt_spec_labels,
                        filename="TEMP_Cu_fcc_v2.dat" )

end

test_Cu_fcc()
