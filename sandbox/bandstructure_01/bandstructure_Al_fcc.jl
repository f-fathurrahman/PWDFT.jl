using Printf
using LinearAlgebra
using Random
using PWDFT

include("dump_bandstructure.jl")
include("gen_kpath.jl")

function test_Al_fcc()

    # Atoms
    atoms = Atoms( xyz_string_frac=
        """
        1

        Al  0.0  0.0  0.0
        """, in_bohr=true,
        LatVecs = gen_lattice_fcc(4.0495*ANG2BOHR) )

    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/Al-q3.gth"]
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc,
                       meshk=[8,8,8], extra_states=4 )
    println(Ham)

    # Solve the KS problem
    KS_solve_SCF!( Ham, betamix=0.2, mix_method="rpulay", use_smearing=true, kT=0.001 )


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
                               Nstates_empty=4 )

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
                        filename="TEMP_Al_fcc_v2.dat" )

end

test_Al_fcc()
