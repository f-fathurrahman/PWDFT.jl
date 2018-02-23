using PWDFT

function test_main()
    #
    atoms = init_atoms_xyz("H.xyz")
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
    #
    update!(Ham, rhoe)
    #
    g = calc_grad(Ham, psi)
    s = sum(g)
    @printf("sum g = (%18.10f,%18.10fim)\n", s.re, s.im)
end

test_main()
