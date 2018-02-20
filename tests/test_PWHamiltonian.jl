using PWDFT

include("PWHamiltonian.jl")
include("ortho_gram_schmidt.jl")
include("calc_rhoe.jl")

function test_main()
    #
    atoms = init_atoms_xyz("H.xyz")
    println(atoms)
    #
    LatVecs = 16.0*diagm(ones(3))
    pw = PWGrid(30.0*0.5, LatVecs)
    println(pw)
    #
    Ham = PWHamiltonian( pw, atoms )
    #
    Ngwx = Ham.pw.gvecw.Ngwx
    Nstates = 1
    Focc = [1.0]
    #
    srand(1234)
    psi = rand(Ngwx,Nstates) + im*rand(Ngwx,Nstates)
    psi = ortho_gram_schmidt(psi)
    #
    rhoe = calc_rhoe( pw, Focc, psi )
    @printf("Integ rhoe = %18.10f\n", sum(rhoe)*pw.Ω/prod(pw.Ns))

    println("\nBefore updating Hamiltonian")
    @time Kpsi = op_K(Ham, psi)
    println(sum(Kpsi))
    #
    @time Vpsi = op_V_loc(Ham, psi)
    println(sum(Vpsi))
    #
    @time Vpsi = op_V_Ps_loc(Ham, psi)
    println(sum(Vpsi))

    update!(Ham, rhoe)

    println("\nAfter updating Hamiltonian")
    @time Kpsi = op_K(Ham, psi)
    println(sum(Kpsi))
    #
    @time Vpsi = op_V_loc(Ham, psi)
    println(sum(Vpsi))
    #
    @time Vpsi = op_V_Ps_loc(Ham, psi)
    println(sum(Vpsi))

end

test_main()
