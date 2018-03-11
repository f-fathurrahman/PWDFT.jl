using PWDFT

function test_main()

    # Atoms
    atoms = init_atoms_xyz("PtO.xyz")
    println(atoms)

    pspfiles = ["../pseudopotentials/pade_gth/Pt-q10.gth",
                "../pseudopotentials/pade_gth/O-q6.gth"]

    ecutwfc = 60.0
    LatVecs = 20.0*diagm(ones(3))
    Ham = PWHamiltonian( atoms, pspfiles, ecutwfc, LatVecs )

    println(Ham.pw)
    for isp = 1:atoms.Nspecies
        println(Ham.pspots[isp])
    end

    Nstates = 4
    Ngwx = Ham.pw.gvecw.Ngwx
    srand(1234)
    psi = rand(Complex128,Ngwx,Nstates)
    psi = ortho_gram_schmidt(psi)

    if Ham.pspotNL.NbetaNL > 0
        betaNL_psi = calc_betaNL_psi( Ham.pspotNL.betaNL, psi )
        E_ps_NL = calc_E_Ps_nloc( Ham, psi )
    end

end

test_main()
