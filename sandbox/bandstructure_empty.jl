using Printf
using LinearAlgebra
using Random
using PWDFT

include("dump_bandstructure.jl")

function test_empty_lattice(lattice::String, band_file::String)
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        1

        H  0.0   0.0   0.0
        """, in_bohr=true)
    if lattice == "sc"
        atoms.LatVecs = gen_lattice_sc(6.0)
    elseif lattice == "fcc"
        atoms.LatVecs = gen_lattice_fcc(6.0)
    elseif lattice == "bcc"
        atoms.LatVecs = gen_lattice_bcc(6.0)
    elseif lattice == "hexagonal"
        atoms.LatVecs = gen_lattice_hexagonal(6.0,10.0)
    elseif lattice == "orthorhombic"
        atoms.LatVecs = gen_lattice_orthorhombic(6.0,10.0,8.0)
    elseif lattice == "monoclinic"
        atoms.LatVecs = gen_lattice_monoclinic(6.0,7.0,6.5,80.0)
    elseif lattice == "tetragonal"
        atoms.LatVecs = gen_lattice_tetragonal_P(6.0,7.0)
    elseif lattice == "rhombohedral"
        atoms.LatVecs = gen_lattice_rhombohedral(6.0,80.0)
    else
        @printf("lattice is not known: %s\n", lattice)
    end
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)

    if lattice == "sc"
        kpoints, kpt_spec, kpt_spec_labels = kpath_from_file(atoms, "KPATH_SC_60")
    elseif lattice == "fcc"
        kpoints, kpt_spec, kpt_spec_labels = kpath_from_file(atoms, "KPATH_FCC_60")
    elseif lattice == "bcc"
        kpoints, kpt_spec, kpt_spec_labels = kpath_from_file(atoms, "KPATH_BCC_60")
    elseif lattice == "hexagonal"
        kpoints, kpt_spec, kpt_spec_labels = kpath_from_file(atoms, "KPATH_HEX_60")
    elseif lattice == "orthorhombic"
        kpoints, kpt_spec, kpt_spec_labels = kpath_from_file(atoms, "KPATH_ORTHO_60")
    elseif lattice == "monoclinic"
        kpoints, kpt_spec, kpt_spec_labels = kpath_from_file(atoms, "KPATH_MONOCLINIC_60")
    elseif lattice == "tetragonal"
        kpoints, kpt_spec, kpt_spec_labels = kpath_from_file(atoms, "KPATH_TETRAGONAL_60")
    elseif lattice == "rhombohedral"
        kpoints, kpt_spec, kpt_spec_labels = kpath_from_file(atoms, "KPATH_RHOMBOHEDRAL_T1_SEG_1_60")
    else
        @printf("lattice is not known: %s\n", lattice)
    end
    println(kpoints)

    ecutwfc = 15.0
    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth"]
    Ham = free_electron_Hamiltonian(
            atoms, pspfiles, ecutwfc, kpoints=kpoints, extra_states=20,
          )

    pw = Ham.pw
    Nkpt = kpoints.Nkpt
    Nspin = Ham.electrons.Nspin
    Nstates = Ham.electrons.Nstates

    Nkspin = Nkpt*Nspin
    Ngw = pw.gvecw.Ngw

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
        
        #evals[:,ikspin], psiks[ikspin] =
        #diag_davidson( Ham, psiks[ikspin], verbose_last=true )

        #evals[:,ikspin], psiks[ikspin] =
        #diag_Emin_PCG( Ham, psiks[ikspin], verbose_last=true )

    end
    end

    dump_bandstructure( evals, kpoints.k, kpt_spec, kpt_spec_labels, filename=band_file )

end

#test_empty_lattice("sc", "TEMP_empty_lattice_sc.dat")
test_empty_lattice("fcc", "TEMP_empty_lattice_fcc.dat")

#test_empty_lattice("bcc", "TEMP_empty_lattice_bcc.dat")

#test_empty_lattice("hexagonal", "TEMP_empty_lattice_hex.dat")

#test_empty_lattice("orthorhombic", "TEMP_empty_lattice_ortho.dat")

#test_empty_lattice("monoclinic", "TEMP_empty_lattice_monoclinic.dat")

#test_empty_lattice("tetragonal", "TEMP_empty_lattice_tetragonal.dat")

#test_empty_lattice("rhombohedral", "TEMP_empty_lattice_rhombohedral.dat")
