using Random
using Printf
using PWDFT

function test_Nkpt_1()
    #
    #atoms = init_atoms_xyz("../structures/H.xyz")
    atoms = init_atoms_xyz("../structures/H2.xyz")
    atoms.LatVecs = gen_lattice_sc(16.0)
    println(atoms)
    #
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, ecutwfc )
    #
    pw = Ham.pw
    Ngwx = pw.gvecw.Ngwx
    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc
    Nkpt = pw.gvecw.kpoints.Nkpt
    Ngw = pw.gvecw.Ngw
    #
    Random.seed!(1234)
    psik = Array{Array{ComplexF64,2},1}(undef,Nkpt)
    for ik = 1:Nkpt
        psi = rand(Ngw[ik],Nstates) + im*rand(Ngw[ik],Nstates)
        psik[ik] = ortho_gram_schmidt(psi)
    end
    #
    rhoe = calc_rhoe( pw, Focc, psik )
    @printf("Integ rhoe = %18.10f\n", sum(rhoe)*pw.CellVolume/prod(pw.Ns))

    println("\nTesting op_K")
    #for ik = 1:Nkpt
    #    Ham.ik = ik  # don't forget set current kpoints index
        @time Kpsi = op_K(Ham, psik[1])
        s = sum(Kpsi)
    #    @printf("ik=%d sum Kpsi = (%18.10f,%18.10fim)\n", ik, s.re, s.im)
    #end

    #
    #for ik = 1:Nkpt
    #    Ham.ik = ik
        @time Vpsi1 = op_V_loc(Ham, psik[1])
    #    s = sum(Vpsi1)
    #    @printf("sum Vpsi1 = (%18.10f,%18.10fim)\n", s.re, s.im)
    #end

    """
    #
    @time Vpsi2 = op_V_Ps_loc(Ham, psi)
    s = sum(Vpsi2)
    @printf("sum Vpsi2 = (%18.10f,%18.10fim)\n", s.re, s.im)
    #
    @time Hpsi = op_H( Ham, psi )
    s = sum(Hpsi - Kpsi - Vpsi1)
    @printf("diff op_H = (%18.10f,%18.10fim)\n", s.re, s.im)

    update!(Ham, rhoe)

    println("\nAfter updating Hamiltonian")
    @time Kpsi = op_K(Ham, psi)
    s = sum(Kpsi)
    @printf("sum Kpsi = (%18.10f,%18.10fim)\n", s.re, s.im)
    #
    @time Vpsi1 = op_V_loc(Ham, psi)
    s = sum(Vpsi1)
    @printf("sum Vpsi1 = (%18.10f,%18.10fim)\n", s.re, s.im)
    #
    @time Vpsi2 = op_V_Ps_loc(Ham, psi)
    s = sum(Vpsi2)
    @printf("sum Vpsi2 = (%18.10f,%18.10fim)\n", s.re, s.im)
    #
    @time Hpsi = op_H( Ham, psi )
    s = sum(Hpsi - Kpsi - Vpsi1)
    @printf("diff op_H = (%18.10f,%18.10fim)\n", s.re, s.im)

    Energies = calc_energies(Ham, psi)
    println("\nCalculated energies")
    println(Energies)

    println("\nOld energies of Hamiltonian:")
    println(Ham.energies)

    Ham.energies = Energies
    println("\nUpdated energies of Hamiltonian:")
    println(Ham.energies)
    """

end

function test_Si()
    #
    atoms = init_atoms_xyz_string(
        """
        2

        Si  0.0   0.0   0.0
        Si  0.25  0.25  0.25
        """, in_bohr=true)
    atoms.LatVecs = gen_lattice_fcc(5.431*ANG2BOHR)
    atoms.positions = atoms.LatVecs*atoms.positions
    println(atoms)
    #
    pspfiles = ["../pseudopotentials/pade_gth/Si-q4.gth"]
    ecutwfc = 25.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3] )
    #
    pw = Ham.pw
    Ngwx = pw.gvecw.Ngwx
    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc
    Nkpt = pw.gvecw.kpoints.Nkpt
    Ngw = pw.gvecw.Ngw
    #
    Random.seed!(1234)
    psik = Array{Array{ComplexF64,2},1}(undef,Nkpt)
    for ik = 1:Nkpt
        psi = rand(Ngw[ik],Nstates) + im*rand(Ngw[ik],Nstates)
        psik[ik] = ortho_gram_schmidt(psi)
    end
    #
    rhoe = calc_rhoe( pw, Focc, psik )
    @printf("Integ rhoe = %18.10f\n", sum(rhoe)*pw.CellVolume/prod(pw.Ns))

    println("\nTesting op_K")
    for ik = 1:Nkpt
        Ham.ik = ik  # don't forget set current kpoints index
        @time Kpsi = op_K(Ham, psik[ik])
        s = sum(Kpsi)
        @printf("ik=%d sum Kpsi = (%18.10f,%18.10fim)\n", ik, s.re, s.im)
    end

    println("\nTesting op_V_Ps_loc")
    for ik = 1:Nkpt
        Ham.ik = ik
        @time Vpsi2 = op_V_Ps_loc(Ham, psik[ik])
        s = sum(Vpsi2)
        @printf("sum Vpsi2 = (%18.10f,%18.10fim)\n", s.re, s.im)
    end

    println("\nTesting op_V_loc")
    for ik = 1:Nkpt
        Ham.ik = ik
        @time Vpsi1 = op_V_loc(Ham, psik[ik])
        s = sum(Vpsi1)
        @printf("sum Vpsi1 = (%18.10f,%18.10fim)\n", s.re, s.im)
    end

    println("\nTesting op_H")
    for ik = 1:Nkpt
        Ham.ik = ik
        @time Hpsi = op_H(Ham, psik[ik])
        s = sum(Hpsi)
        @printf("sum Hpsi = (%18.10f,%18.10fim)\n", s.re, s.im)
    end

end


function test_Ni_fcc(;Nspin=1)

    @assert( Nspin <= 2 )

    #
    atoms = init_atoms_xyz_string(
        """
        1

        Ni  0.0   0.0   0.0
        """)
    atoms.LatVecs = gen_lattice_fcc(6.65914911201)
    println(atoms)
    #
    pspfiles = ["../pseudopotentials/pade_gth/Ni-q10.gth"]
    ecutwfc = 25.0
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], Nspin=Nspin )

    #
    pw = Ham.pw
    Ngwx = pw.gvecw.Ngwx
    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc
    Nkpt = pw.gvecw.kpoints.Nkpt
    Ngw = pw.gvecw.Ngw
    Nkspin = Nkpt*Nspin
    Npoints = prod(pw.Ns)
    Ng = pw.gvec.Ng
    G = pw.gvec.G
    G2 = pw.gvec.G2
    idx_g2r = pw.gvec.idx_g2r
    CellVolume = pw.CellVolume

    #
    Random.seed!(1234)
    psiks = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        ikspin = ik + (ispin - 1)*Nkpt
        psi = rand(Ngw[ik],Nstates) + im*rand(Ngw[ik],Nstates)
        psiks[ikspin] = ortho_gram_schmidt(psi)
    end
    end

    Rhoe = zeros(Npoints,Nspin)
    dVol = pw.CellVolume/prod(pw.Ns)
    for ispin = 1:Nspin
        idxset = (Nkpt*(ispin-1)+1):(Nkpt*ispin)
        Rhoe[:,ispin] = calc_rhoe( pw, Focc[:,idxset], psiks[idxset], renormalize=false )
        @printf("spin %d Integ Rhoe = %18.10f\n", ispin, sum(Rhoe[:,ispin])*dVol)
    end

    update!(Ham, Rhoe)

    println("\nTesting op_K")
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        Ham.ik = ik  # don't forget set current kpoints index
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        Kpsi = op_K(Ham, psiks[ikspin])
        s = sum(Kpsi)
        @printf("ispin = %d ik=%d sum Kpsi = (%18.10f,%18.10fim)\n", ispin, ik, s.re, s.im)
    end
    end

    println("\nTesting op_V_Ps_loc")
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        Vpsi2 = op_V_Ps_loc(Ham, psiks[ikspin])
        s = sum(Vpsi2)
        @printf("ispin = %d ik=%d sum Vpsi2 = (%18.10f,%18.10fim)\n", ispin, ik, s.re, s.im)
    end # ik
    end # ispin

    println("\nTesting op_V_loc")
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        Vpsi1 = op_V_loc(Ham, psiks[ikspin])
        s = sum(Vpsi1)
        @printf("ispin = %d ik=%d sum Vpsi1 = (%18.10f,%18.10fim)\n", ispin, ik, s.re, s.im)
    end # ik
    end # ispin

    println("\nTesting op_H")
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt        
        Hpsi = op_H(Ham, psiks[ikspin])
        s = sum(Hpsi)
        @printf("ispin = %d ik=%d sum Hpsi = (%18.10f,%18.10fim)\n", ispin, ik, s.re, s.im)        
    end # ik
    end # ispin

    Energies = calc_energies( Ham, psiks )
    println(Energies)

    if Nspin == 2
        Rhoe_total = Rhoe[:,1] + Rhoe[:,2]
    else
        Rhoe_total = Rhoe[:,1]
    end
    rhoG = R_to_G( pw, Rhoe_total ) / Npoints  # XXX normalize by Npoints
    
    phiG = Poisson_solve( pw, Rhoe_total)
    #phiG = R_to_G( pw, Ham.potentials.Hartree )
    EhartreeG = 0.0
    for ig = 2:Ng
        ip = idx_g2r[ig]
        EhartreeG = EhartreeG + 0.5*real(phiG[ip]*conj(rhoG[ip]))*CellVolume/Npoints
    end
    println("EhartreeG = ", EhartreeG)
    println("diff Ehartree = ", abs(EhartreeG - Energies.Hartree))

    # assumption: use VWN
    epsxcG = R_to_G( pw, calc_epsxc_VWN( Rhoe ) )
    ExcG = 0.0
    for ig = 1:Ng  # ig=1 is included
        ip = idx_g2r[ig]
        ExcG = ExcG + real(epsxcG[ip]*conj(rhoG[ip]))*CellVolume/Npoints
    end
    println("ExcG = ", ExcG)
    println("diff Exc = ", abs(ExcG - Energies.XC))

    V_ps_loc_G = R_to_G( pw, Ham.potentials.Ps_loc )
    E_ps_loc_G = 0.0
    for ig = 1:Ng  # ig=1 is included
        ip = idx_g2r[ig]
        E_ps_loc_G = E_ps_loc_G + real(V_ps_loc_G[ip]*conj(rhoG[ip]))*CellVolume/Npoints
    end
    println("E_ps_log_G = ", E_ps_loc_G)
    println("diff E_ps_loc = ", abs(E_ps_loc_G - Energies.Ps_loc))

end

test_Si()

test_Ni_fcc(Nspin=1)
test_Ni_fcc(Nspin=2)

test_Nkpt_1()
