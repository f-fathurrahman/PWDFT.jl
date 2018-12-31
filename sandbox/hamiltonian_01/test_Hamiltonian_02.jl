using Random
using Printf
using PWDFT

function init_random_psik( Ham::Hamiltonian )
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Ngw = Ham.pw.gvecw.Ngw
    Nstates = Ham.electrons.Nstates
    return init_random_psik( Ngw, Nstates, Nkpt )
end

function init_random_psik( Ngw::Array{Int64,1}, Nstates::Int64, Nkpt::Int64 )
    @assert( length(Ngw) == Nkpt )
    psik = Array{Array{ComplexF64,2},1}(undef,Nkpt)
    for ik = 1:Nkpt
        psik[ik] = rand(ComplexF64,Ngw[ik],Nstates)
        ortho_gram_schmidt!(psik[ik])
    end
    return psik
end

function init_random_psi( Ham::Hamiltonian )
    @assert( Ham.electrons.Nspin == 1 )
    @assert( Ham.pw.gvecw.kpoints.Nkpt == 1 )
    return init_random_psi( Ham.pw.gvecw.Ngw[1], Ham.electrons.Nstates )
end

function init_random_psi( Ngw_ik::Int64, Nstates::Int64 )
    return ortho_gram_schmidt( rand(ComplexF64,Ngw_ik,Nstates) )
end

function test_Nkpt_1()
    
    atoms = init_atoms_xyz("../structures/O2.xyz")
    #atoms = init_atoms_xyz("../structures/H2.xyz")
    atoms.LatVecs = gen_lattice_sc(16.0)
    
    ecutwfc = 15.0
    Ham = Hamiltonian( atoms, ecutwfc )
    
    # Shortcuts
    pw = Ham.pw
    Focc = Ham.electrons.Focc
    CellVolume = pw.CellVolume
    Ns = pw.Ns
    
    Random.seed!(4321)
    psi = init_random_psi( Ham )

    rhoe = calc_rhoe( Ham, psi )
    @printf("Integ rhoe = %18.10f\n", sum(rhoe)*CellVolume/prod(Ns))

    println("\nTesting * operator")

    Hpsi_1 = op_H( Ham, psi )
    Hpsi_2 = Ham * psi

    println("sum Hpsi_1 = ", sum(Hpsi_1))
    println("sum Hpsi_2 = ", sum(Hpsi_2))

    psiH_1 = adjoint( op_H(Ham, psi) )
    psiH_2 = psi'*Ham

    println("sum psiH_1 = ", sum(psiH_1))
    println("sum psiH_2 = ", sum(psiH_2))

    println("\nTesting adjoint * operator")

    Hsub_1 = psiH_2 * psi
    println("sum Hsub_1 = ", sum(Hsub_1))

    Hsub_2 = psi'*Ham * psi
    println("sum Hsub_2 = ", sum(Hsub_2))

    Hsub_3 = psi' * Hpsi_1
    println("sum Hsub_3 = ", sum(Hsub_3))

    psi1 = init_random_psi(Ham)
    psi2 = init_random_psi(Ham)

    psi2T = psi2'

    println("\nThis gives small imaginary component")
    println("diff = ", sum(psi2'*op_H(Ham,psi1) - psi2'*Ham * psi1))

    println("\nThis should be exactly zero")
    println("diff = ", sum(psi2'*op_H(Ham,psi1) - psi2'*(Ham * psi1)))
    
    println("\nThis should be exactly zero")
    println("diff = ", sum(adjoint(op_H(Ham,psi2T'))*psi1 - psi2'*Ham * psi1))

    println("\nThis gives small imaginary component")
    println("diff = ", sum(adjoint(op_H(Ham,psi2T'))*psi1 - psi2'*(Ham * psi1)))

end

test_Nkpt_1()
