using PWDFT

function test_main()
    #
    atoms = init_atoms_xyz("../structures/H.xyz")
    println(atoms)
    #
    LatVecs = 16.0*diagm(ones(3))
    pw = PWGrid(30.0, LatVecs)
    println(pw)
    #
    Ham = PWHamiltonian( pw, atoms )
    #
    Ngwx = Ham.pw.gvecw.Ngwx
    Nstates = 1
    Focc = [1.0]
    Ham.focc = Focc
    #
    srand(1234)
    psi = rand(Ngwx,Nstates) + im*rand(Ngwx,Nstates)
    psi = ortho_gram_schmidt(psi)
    #
    rhoe = calc_rhoe( pw, Focc, psi )
    @printf("Integ rhoe = %18.10f\n", sum(rhoe)*pw.Ω/prod(pw.Ns))

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
end

test_main()
