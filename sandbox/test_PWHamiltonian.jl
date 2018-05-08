using PWDFT

function test_Nkpt_1()
    #
    #atoms = init_atoms_xyz("../structures/H.xyz")
    atoms = init_atoms_xyz("../structures/H2.xyz")
    println(atoms)
    #
    ecutwfc = 15.0
    Ham = PWHamiltonian( atoms, ecutwfc, verbose=true )
    #
    pw = Ham.pw
    Ngwx = pw.gvecw.Ngwx
    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc
    Nkpt = pw.gvecw.kpoints.Nkpt
    Ngw = pw.gvecw.Ngw
    #
    srand(1234)
    psik = Array{Array{Complex128,2},1}(Nkpt)
    for ik = 1:Nkpt
        psi = rand(Ngw[ik],Nstates) + im*rand(Ngw[ik],Nstates)
        psik[ik] = ortho_gram_schmidt(psi)
    end
    #
    rhoe = calc_rhoe( pw, Focc, psik )
    @printf("Integ rhoe = %18.10f\n", sum(rhoe)*pw.Ω/prod(pw.Ns))

    println("\nTesting op_K")
    for ik = 1:Nkpt
        Ham.ik = ik  # don't forget set current kpoints index
        @time Kpsi = op_K(Ham, psik[ik])
        s = sum(Kpsi)
        @printf("ik=%d sum Kpsi = (%18.10f,%18.10fim)\n", ik, s.re, s.im)
    end

    #
    for ik = 1:Nkpt
        Ham.ik = ik
        @time Vpsi1 = op_V_loc(Ham, psi)
        s = sum(Vpsi1)
        @printf("sum Vpsi1 = (%18.10f,%18.10fim)\n", s.re, s.im)
    end

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
    Ham = PWHamiltonian( atoms, pspfiles, ecutwfc, meshk=[3,3,3], verbose=true )
    #
    pw = Ham.pw
    Ngwx = pw.gvecw.Ngwx
    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc
    Nkpt = pw.gvecw.kpoints.Nkpt
    Ngw = pw.gvecw.Ngw
    #
    srand(1234)
    psik = Array{Array{Complex128,2},1}(Nkpt)
    for ik = 1:Nkpt
        psi = rand(Ngw[ik],Nstates) + im*rand(Ngw[ik],Nstates)
        psik[ik] = ortho_gram_schmidt(psi)
    end
    #
    rhoe = calc_rhoe( pw, Focc, psik )
    @printf("Integ rhoe = %18.10f\n", sum(rhoe)*pw.Ω/prod(pw.Ns))

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

test_Si()

#test_Nkpt_1()
