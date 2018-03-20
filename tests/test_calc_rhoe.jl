using PWDFT

function test_main()
    ecutwfc_Ry = 30.0
    LatVecs = 16.0*diagm(ones(3))
    pw = PWGrid( 0.5*ecutwfc_Ry, LatVecs )

    Nstates = 4
    Focc = 2.0*ones(Nstates)

    Ngwx = pw.gvecw.Ngwx
    psi = rand( Complex128, Ngwx, Nstates )
    psi = ortho_gram_schmidt(psi)  # orthogonalize in G-space

    rhoe = calc_rhoe( pw, Focc, psi )
    dVol = pw.Ω/prod(pw.Ns)
    @printf("Integrated rhoe = %18.10f\n", sum(rhoe)*dVol)

    rhoeG_full = R_to_G( pw, rhoe )
    println("sum rhoeG_full = ", sum(rhoeG_full))

    Ng = pw.gvec.Ng
    rhoeG = zeros( Complex128, Ng )
    rhoeG = rhoeG_full[pw.gvec.idx_g2r]
    println("sum rhoeG = ", sum(rhoeG))
end

test_main()

