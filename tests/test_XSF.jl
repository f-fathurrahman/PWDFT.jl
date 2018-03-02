using PWDFT

function test_main()
    LatVecs = 16.0*diagm( ones(3) )
    ecutRy = 30.0
    pw = PWGrid( 0.5*ecutRy, LatVecs )
    println( pw )

    Natoms = 1
    atpos = zeros( 3, Natoms )
    atpos[:,1] = [ 8.0, 8.0, 7.0 ]

    Npoints = prod( pw.Ns )
    idx = pw.gvecw.idx_gw2r
    gwave = pw.gvec.G[:,idx]
    Ngwx = pw.gvecw.Ngwx

    l = 2
    ia = 1
    for m = -l:l
        ctmp = zeros( Complex128, Npoints )
        psi = zeros( Complex128, Ngwx )
        for ig = 1:Ngwx
            g = gwave[:,ig]
            Gm = norm(g)
            GX = atpos[1,ia]*g[1] + atpos[2,ia]*g[2] + atpos[3,ia]*g[3]
            Sf = cos(GX) - im*sin(GX)
            psi[ig] = (1.0*im)^l * Ylm_real(l,m,g) * exp( -0.5*Gm^2 ) * Sf
        end
        ctmp[idx] = psi
        psiR_real = real( G_to_R( pw, ctmp ) ) * Npoints
        filxsf = "psi_l_" * string(l) * "_" * string(m) * ".xsf"
        write_xsf( filxsf, LatVecs/ANG2BOHR, atpos/ANG2BOHR, molecule=false )
        write_xsf_data3d_crystal(  filxsf, pw.Ns, LatVecs/ANG2BOHR, psiR_real )
        @printf("File %s is written\n", filxsf)
    end

end

test_main()


