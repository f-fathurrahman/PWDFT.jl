using PWDFT

function test_main()
    pw = PWGrid(60.0, 20.0*diagm(ones(3)))
    println(pw)

    # Atoms
    atoms = init_atoms_xyz("Pt.xyz")
    println(atoms)

    # Structure factor
    strf = calc_strfact( atoms, pw )

    Nspecies = 1
    Pspots = Array{PsPot_GTH}(Nspecies)
    Pspots[1] = PsPot_GTH("../pseudopotentials/pade_gth/Pt-q10.gth")

    @time pspotNL = PsPotNL( pw, atoms, Pspots )

    Nstates = 4
    Ngwx = pw.gvecw.Ngwx
    srand(1234)
    psi = rand(Complex128,Ngwx,Nstates)
    psi = ortho_gram_schmidt(psi)

    betaNL_psi = calc_betaNL_psi( pspotNL.betaNL, psi )

end

test_main()

