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
            atoms, pspfiles, 15.0, kpoints=kpoints, extra_states=3,
            verbose=true
          )

    pw = Ham.pw
    Nkpt = kpoints.Nkpt
    Nspin = Ham.electrons.Nspin
    Nstates = Ham.electrons.Nstates

    Nkspin = Nkpt*Nspin
    Ngw = pw.gvecw.Ngw

    psiks = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    evals = zeros(Float64,Nstates,Nkspin)

    srand(1234)
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
        
        #evals[:,ikspin], psiks[ikspin] =
        #diag_davidson( Ham, psiks[ikspin], verbose_last=true )

        evals[:,ikspin], psiks[ikspin] =
        diag_Emin_PCG( Ham, psiks[ikspin], verbose_last=true )

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

function save_potential(Ham::Hamiltonian)
    # TO BE IMPLEMENTED
end


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
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5, meshk=[3,3,3], verbose=true )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

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
                               Nstates_empty=2 )

    Nstates = Ham.electrons.Nstates
    ebands = Ham.electrons.ebands
    Nspin = Ham.electrons.Nspin

    Nkspin = Nkpt*Nspin

    # Also, new PsPotNL
    Ham.pspotNL = PsPotNL( pw, atoms, Ham.pspots, kpoints, check_norm=false )


    psiks = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    evals = zeros(Float64,Nstates,Nkspin)

    srand(1234)
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        ikspin = ik + (ispin - 1)*Nkpt
        psiks[ikspin] = ortho_gram_schmidt(rand(ComplexF64,Ngw[ik],Nstates))
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
        diag_davidson( Ham, psiks[ikspin], verbose=true )

    end
    end

    dump_bandstructure( evals, kpoints.k, kpt_spec, kpt_spec_labels, filename="TEMP_Si_fcc.dat" )


end

#test_Si_fcc()




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
    ecutwfc_Ry = 30.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc_Ry*0.5,
                         meshk=[8,8,8], verbose=true, extra_states=4 )

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( atoms )

    #
    # Solve the KS problem
    #
    KS_solve_SCF_smearing!( Ham, β=0.2, mix_method="anderson" )


    #
    # Band structure calculation
    #


    kpoints, kpt_spec, kpt_spec_labels = kpath_from_file(atoms, "KPATH_FCC_60")

    # New pw
    pw = PWGrid(15.0, atoms.LatVecs, kpoints=kpoints)
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


    psiks = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    evals = zeros(Float64,Nstates,Nkspin)

    srand(1234)
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        ikspin = ik + (ispin - 1)*Nkpt
        psiks[ikspin] = ortho_gram_schmidt(rand(ComplexF64,Ngw[ik],Nstates))
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

    dump_bandstructure( evals, kpoints.k, kpt_spec, kpt_spec_labels, filename="TEMP_Cu_fcc.dat" )

end

#test_Cu_fcc()
