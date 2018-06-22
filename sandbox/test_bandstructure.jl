using Printf
using LinearAlgebra
using Random
using PWDFT

include("dump_bandstructure.jl")


function test_bandstructure_fcc()
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        1

        H  0.0   0.0   0.0
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(6.0)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)
    #
    kpoints, kpt_spec, kpt_spec_labels = kpath_from_file(atoms, "KPATH_FCC_60")
    println(kpoints)
    #
    pw = PWGrid(15.0, atoms.LatVecs, kpoints=kpoints)
    #
    # Prepare PWHamiltonian
    #
    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth"]
    #
    Nspecies = atoms.Nspecies
    if Nspecies != size(pspfiles)[1]
        @printf("ERROR length of pspfiles is not equal to %d\n", Nspecies)
        exit()
    end

    Npoints = prod(pw.Ns)
    Ω = pw.Ω
    G2 = pw.gvec.G2
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r

    strf = calc_strfact( atoms, pw )

    #
    # Initialize pseudopotentials and local potentials
    #
    Vg = zeros(ComplexF64, Npoints)
    V_Ps_loc = zeros(Float64, Npoints)

    Pspots = Array{PsPot_GTH}(undef,Nspecies)

    for isp = 1:Nspecies
        Pspots[isp] = PsPot_GTH( pspfiles[isp] )
        psp = Pspots[isp]
        println(psp)
        for ig = 1:Ng
            ip = idx_g2r[ig]
            #Vg[ip] = strf[ig,isp] * eval_Vloc_G( psp, G2[ig] ) / Ω
            Vg[ip] = 0.0
        end
        #
        V_Ps_loc[:] = V_Ps_loc[:] + real( G_to_R(pw, Vg) ) * Npoints
    end

    Nspin = 1

    # other potential terms are set to zero
    V_Hartree = zeros( Float64, Npoints )
    V_XC = zeros( Float64, Npoints, Nspin )
    potentials = Potentials( V_Ps_loc, V_Hartree, V_XC )
    #
    energies = Energies()
    #
    rhoe = zeros( Float64, Npoints, Nspin )

    #electrons = Electrons( atoms, Pspots, Nspin=Nspin, Nkpt=kpoints.Nkpt,
    #                       Nstates_empty=0 )
    #println(electrons)
    Nkpt = kpoints.Nkpt
    #
    electrons = Electrons(
        8, #Nelectrons::Float64
        4, #Nstates::Int64
        4, #Nstates_occ::Int64
        2*ones(4,Nkpt), #Focc::Array{Float64,2}
        zeros(4,Nkpt), #ebands::Array{Float64,2}
        1 #Nspin::Int64
    )

    # NL pseudopotentials
    pspotNL = PsPotNL()

    atoms.Zvals = get_Zvals( Pspots )

    ik = 1
    ispin = 1
    xcfunc = "VWN"
    Ham = PWHamiltonian( pw, potentials, energies, rhoe,
                         electrons, atoms, Pspots, pspotNL, xcfunc, ik, ispin )

    Nkspin = Nkpt*Nspin
    Ngw = pw.gvecw.Ngw
    Nstates = electrons.Nstates

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
        @printf("\nispin = %d, ik = %d, ikspin=%d\n", ispin, ik, ikspin)
        @printf("kpts = [%f,%f,%f]\n", k[1,ik], k[2,ik], k[3,ik])
        #
        #evals[:,ikspin], psiks[ikspin] =
        #diag_davidson( Ham, psiks[ikspin], verbose_last=true, NiterMax=300 )
        #
        #if ikspin == 55 # problematic kpoints
        #    evals[:,ikspin], psiks[ikspin] =
        #    diag_Emin_PCG( Ham, psiks[ikspin], verbose=true, NiterMax=300 )
        #else
            evals[:,ikspin], psiks[ikspin] =
            diag_lobpcg( Ham, psiks[ikspin], verbose_last=true )
        #end
        #
    end
    end

    dump_bandstructure( evals, kpoints.k, kpt_spec, kpt_spec_labels, filename="TEMP_band_free_fcc.dat" )

end


function test_bandstructure_sc()
    # Atoms
    atoms = init_atoms_xyz_string(
        """
        1

        H  0.0   0.0   0.0
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_sc(5.0)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)
    #
    kpoints, kpt_spec, kpt_spec_labels = kpath_from_file(atoms, "KPATH_SC_60")
    println(kpoints)
    #
    pw = PWGrid(15.0, atoms.LatVecs, kpoints=kpoints)
    println(pw)
    #
    # Prepare PWHamiltonian
    #
    pspfiles = ["../pseudopotentials/pade_gth/H-q1.gth"]
    #
    Nspecies = atoms.Nspecies
    if Nspecies != size(pspfiles)[1]
        @printf("ERROR length of pspfiles is not equal to %d\n", Nspecies)
        exit()
    end

    Npoints = prod(pw.Ns)
    Ω = pw.Ω
    G2 = pw.gvec.G2
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r

    strf = calc_strfact( atoms, pw )

    #
    # Initialize pseudopotentials and local potentials
    #
    Vg = zeros(ComplexF64, Npoints)
    V_Ps_loc = zeros(Float64, Npoints)

    Pspots = Array{PsPot_GTH}(undef,Nspecies)

    for isp = 1:Nspecies
        Pspots[isp] = PsPot_GTH( pspfiles[isp] )
        psp = Pspots[isp]
        println(psp)
        for ig = 1:Ng
            ip = idx_g2r[ig]
            #Vg[ip] = strf[ig,isp] * eval_Vloc_G( psp, G2[ig] ) / Ω
            Vg[ip] = 0.0
        end
        #
        V_Ps_loc[:] = V_Ps_loc[:] + real( G_to_R(pw, Vg) ) * Npoints
    end

    Nspin = 1

    # other potential terms are set to zero
    V_Hartree = zeros( Float64, Npoints )
    V_XC = zeros( Float64, Npoints, Nspin )
    potentials = Potentials( V_Ps_loc, V_Hartree, V_XC )
    #
    energies = Energies()
    #
    rhoe = zeros( Float64, Npoints, Nspin )

    #electrons = Electrons( atoms, Pspots, Nspin=Nspin, Nkpt=kpoints.Nkpt,
    #                       Nstates_empty=0 )
    #println(electrons)
    Nkpt = kpoints.Nkpt
    #
    electrons = Electrons(
        8, #Nelectrons::Float64
        4, #Nstates::Int64
        4, #Nstates_occ::Int64
        2*ones(4,Nkpt), #Focc::Array{Float64,2}
        zeros(4,Nkpt), #ebands::Array{Float64,2}
        1 #Nspin::Int64
    )

    # NL pseudopotentials
    pspotNL = PsPotNL()

    atoms.Zvals = get_Zvals( Pspots )

    ik = 1
    ispin = 1
    xcfunc = "VWN"
    Ham = PWHamiltonian( pw, potentials, energies, rhoe,
                         electrons, atoms, Pspots, pspotNL, xcfunc, ik, ispin )

    Nkspin = Nkpt*Nspin
    Ngw = pw.gvecw.Ngw
    Nstates = electrons.Nstates

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
        @printf("\nispin = %d, ik = %d, ikspin=%d\n", ispin, ik, ikspin)
        @printf("kpts = [%f,%f,%f]\n", k[1,ik], k[2,ik], k[3,ik])
        #
        evals[:,ikspin], psiks[ikspin] =
        diag_lobpcg( Ham, psiks[ikspin], verbose_last=true )
    end
    end

    dump_bandstructure( evals, kpoints.k, kpt_spec, kpt_spec_labels, filename="TEMP_band_free_sc.dat" )

end

#test_bandstructure_fcc()
test_bandstructure_sc()
