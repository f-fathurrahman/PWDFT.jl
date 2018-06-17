push!(LOAD_PATH,"../src")

using Printf
using Random
using LinearAlgebra
using PWDFT

include("../src/calc_grad_v2.jl")

function test_main()
    #
    atoms = init_atoms_xyz("../structures/CuSO4.xyz")
    atoms.LatVecs = gen_lattice_sc(16.0)
    println(atoms)
    #
    Ham = PWHamiltonian( atoms, 15.0, verbose=true, Nspin=2 )
    #
    pw = Ham.pw
    Ngwx = pw.gvecw.Ngwx
    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc
    #
    srand(1234)
    psi = ortho_gram_schmidt(rand(ComplexF64,Ngwx,Nstates))
    #
    rhoe = calc_rhoe( 1, pw, Focc, psi )
    @printf("Integ rhoe = %18.10f\n", sum(rhoe)*pw.Ω/prod(pw.Ns))
    #
    update!(Ham, rhoe)
    #
    Nspin = Ham.electrons.Nspin
    for ispin = 1:Nspin
        Ham.ispin = ispin
        # ik is default to 1
        g = calc_grad(Ham, psi)
        s = sum(g)
        @printf("sum g = (%18.10f,%18.10fim)\n", s.re, s.im)
        #
        g = calc_grad_v2(Ham, psi)
        s = sum(g)
        @printf("sum g = (%18.10f,%18.10fim)\n", s.re, s.im)
    end
end

test_main()
