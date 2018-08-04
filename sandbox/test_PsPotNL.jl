using Random
using Printf
using PWDFT

function test_main()

    # Atoms
    atoms = init_atoms_xyz("../structures/PtO.xyz")
    println(atoms)

    pspfiles = ["../pseudopotentials/pade_gth/Pt-q10.gth",
                "../pseudopotentials/pade_gth/O-q6.gth"]

    ecutwfc = 60.0
    LatVecs = gen_lattice_sc(20.0)
    Ham = Hamiltonian( atoms, pspfiles, ecutwfc, verbose=true, extra_states=2 )

    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin
    Nkspin = Nkpt*Nspin
    Nstates = Ham.electrons.Nstates
    Ngw = Ham.pw.gvecw.Ngw
    Random.seed!(1234)
    psiks = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    #
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        ikspin = ik + (ispin-1)*Nkpt
        psiks[ikspin] = ortho_gram_schmidt(rand(ComplexF64,Ngw[ik],Nstates))
    end
    end

    if Ham.pspotNL.NbetaNL > 0
        E_ps_NL = calc_E_Ps_nloc( Ham, psiks )
        println("E ps NL = ", E_ps_NL)
    end
    

end

test_main()
