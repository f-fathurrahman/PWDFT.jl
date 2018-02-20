using PWDFT

include("PWHamiltonian.jl")
include("ortho_gram_schmidt.jl")

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
    Nstates = 4
    #
    srand(1234)
    psi = rand(Ngwx,Nstates) + im*rand(Ngwx,Nstates)
    psi = ortho_gram_schmidt(psi)
    #
    @time Kpsi = op_K(Ham, psi)
    println(sum(Kpsi))
    #
    @time Vpsi = op_V_Ps_loc(Ham, psi)
    println(sum(Vpsi))
end

test_main()
