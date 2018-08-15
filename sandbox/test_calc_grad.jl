using Printf
using Random
using LinearAlgebra
using PWDFT

function test_main()
    #
    atoms = init_atoms_xyz("../structures/H2.xyz")
    atoms.LatVecs = gen_lattice_sc(16.0)
    println(atoms)
    #
    Ham = Hamiltonian( atoms, 15.0, Nspin=2, extra_states=1 )
    #
    pw = Ham.pw
    Npoints = prod(pw.Ns)
    Ngwx = pw.gvecw.Ngwx
    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc
    #
    Nspin = Ham.electrons.Nspin
    Nkpt = 1
    Nkspin = Nkpt*Nspin
    psiks = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    Random.seed!(1234)
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        ikspin = ik + (ispin - 1)*Nkpt
        psiks[ikspin] = ortho_gram_schmidt(rand(ComplexF64,Ngwx,Nstates))
    end
    end
    #
    Rhoe = zeros(Float64,Npoints,Nspin)
    for ispin = 1:Nspin
        idxset = (Nkpt*(ispin-1)+1):(Nkpt*ispin)
        Rhoe[:,ispin] = calc_rhoe( pw, Focc[:,idxset], psiks[idxset] )
    end
    update!(Ham, Rhoe)
    #
    @printf("Integ rhoe = %18.10f\n", sum(Rhoe)*pw.CellVolume/prod(pw.Ns))
    #
    update!(Ham, Rhoe)
    #
    Nspin = Ham.electrons.Nspin
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        ikspin = ik + (ispin-1)*Nkpt
        # ik is default to 1
        g = calc_grad(Ham, psiks[ikspin])
        s = sum(g)
        @printf("%d %d: sum g = (%18.10f,%18.10fim)\n", ispin, ik, s.re, s.im)
    end
    end
end

test_main()
