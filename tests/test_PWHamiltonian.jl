using PWDFT

function test_main()
    #
    atoms = init_atoms_xyz("../structures/H.xyz")
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

"""
    println("\nBefore updating Hamiltonian")
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

test_main()
