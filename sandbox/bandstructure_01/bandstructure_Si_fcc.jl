using Printf
using LinearAlgebra
using Random
using PWDFT

include("dump_bandstructure.jl")

function test_Si_fcc()
    atoms = init_atoms_xyz_string(
        """
        2

        Si  0.0  0.0  0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(10.2631)
    atoms.positions = atoms.LatVecs*atoms.positions
    # Initialize Hamiltonian
    pspfiles = ["../pseudopotentials/pade_gth/Si-q4.gth"]
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, meshk=[3,3,3] )

    KS_solve_Emin_PCG!( Ham, verbose=true )
    
    println("\nTotal energy components")
    println(Ham.energies)

    # Band structure calculation

    kpoints, kpt_spec, kpt_spec_labels = kpath_from_file(atoms, "KPATH_FCC_60")

    # New pw
    pw = PWGrid(15.0, atoms.LatVecs, kpoints=kpoints)
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
        #evals[:,ikspin], psiks[ikspin] =
        #diag_LOBPCG( Ham, psiks[ikspin], verbose_last=true )
        evals[:,ikspin], psiks[ikspin] =
        diag_davidson( Ham, psiks[ikspin], verbose_last=true )

    end
    end

    dump_bandstructure( evals, kpoints.k, kpt_spec, kpt_spec_labels, filename="TEMP_Si_fcc.dat" )

end

test_Si_fcc()
