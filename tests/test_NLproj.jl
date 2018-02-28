using PWDFT

function test_main()
    pw = PWGrid(30.0, 10.0*diagm(ones(3)))
    println(pw)

    # Atoms
    atoms = init_atoms_xyz("Pt.xyz")
    println(atoms)

    # Structure factor
    strf = calc_strfact( atoms, pw )

    zvals = [4.0]
    psp = PsPot_GTH("../pseudopotentials/pade_gth/Pt-q18.gth")
    println(psp)

    Ngwx = pw.gvecw.Ngwx
    idx = pw.gvecw.idx_gw2r
    gwave = pw.gvec.G[:,idx]

    println(size(gwave))

    Nstates = 4
    psi = rand(Ngwx,Nstates) + im*rand(Ngwx,Nstates)
    psi = ortho_gram_schmidt(psi)

    betaNL = zeros( Complex128, Ngwx )

    l = 2
    m = 0
    iprj = 1
    Ω = pw.Ω
    for ig = 1:Ngwx
        g = gwave[:,ig]
        Gm = norm(g)
        betaNL[ig] = Ylm_real( l, m, g ) * eval_proj_G( psp, l, iprj, Gm, Ω )
        @printf("%18.10f %18.10f + %18.10fim\n", Gm, betaNL[ig].re, betaNL[ig].im)
    end
    
    println("rrl = ", psp.rc)
    println("sum betaNL = ", sum(betaNL))

end

test_main()
