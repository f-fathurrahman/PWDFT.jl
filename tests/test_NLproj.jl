using PWDFT

function test_main()
    pw = PWGrid(30.0, 20.0*diagm(ones(3)))
    println(pw)

    # Atoms
    atoms = init_atoms_xyz("Pt.xyz")
    println(atoms)

    # Structure factor
    strf = calc_strfact( atoms, pw )

    Nspecies = 1
    Pspots = Array{PsPot_GTH}(Nspecies)
    Pspots[1] = PsPot_GTH("../pseudopotentials/pade_gth/Pt-q18.gth")

    psp = Pspots[1]
    println(psp)

    Ngwx = pw.gvecw.Ngwx
    idx = pw.gvecw.idx_gw2r
    gwave = pw.gvec.G[:,idx]

    println(size(gwave))

    Natoms = atoms.Natoms

    # 4: indexed from 0:3
    # 0, 1, 2, 3  -> l indexed
    # 1, 2, 3, 4  -> l + 1

    # -3:3
    # -3, -2, -1, 0, 1, 2, 3  -> m
    #  1,  2,  3, 4, 5, 6, 7  -> 4 + m

    # -2, -1, 0, 1, 2  -> m
    #  1,  2, 3, 4, 5  -> 3 + m  lmax = 2 + 1

    prj2beta = Array{Int64}( 3, Natoms, 4, 7 )
    prj2beta[:,:,:,:] = -1   # set to invalid index

    atm2species = atoms.atm2species
    atpos = atoms.positions

    NbetaNL = 0
    for ia = 1:Natoms
        isp = atm2species[ia]
        for l = 0:psp.lmax
            for iprj = 1:psp.Nproj_l[l+1]
                for m = -l:l
                    NbetaNL = NbetaNL + 1
                    prj2beta[iprj,ia,l+1,m+psp.lmax+1] = NbetaNL
                    @printf("NbetaNL, l, m: %3d %3d %3d\n", NbetaNL, l, m)
                end
            end
        end
    end

    betaNL = zeros( Complex128, Ngwx, NbetaNL )

    ibeta = 0
    for ia = 1:Natoms
        isp = atm2species[ia]
        for l = 0:psp.lmax
            for iprj = 1:psp.Nproj_l[l+1]
                for m = -l:l
                    ibeta = ibeta + 1
                    for ig = 1:Ngwx
                        g = gwave[:,ig]
                        Gm = norm(g)
                        GX = atpos[1,ia]*g[1] + atpos[2,ia]*g[2] + atpos[3,ia]*g[3]
                        Sf = cos(GX) - im*sin(GX)
                        betaNL[ig,ibeta] = (1.0*im)^l * Ylm_real(l,m,g) * eval_proj_G( psp, l, iprj, Gm, pw.Ω ) * Sf
                    end
                end
            end
        end
    end

    Npoints = prod(pw.Ns)

    #
    # Check normalization in real space and reciprocal space
    #
    for ibeta = 1:NbetaNL
        norm_G = dot( betaNL[:,ibeta], betaNL[:,ibeta] )
        ctmp = zeros( Complex128, Npoints )
        ctmp[idx] = betaNL[:,ibeta]
        ctmp = G_to_R(pw, ctmp)*Npoints
        data3d = real(ctmp)
        data3d_im = imag(ctmp)
        #
        integ_prj = sum( data3d.^2 ) * pw.Ω / Npoints
        integ_prj_im = sum( data3d_im.^2 ) * pw.Ω / Npoints
        integ_prj_abs = sum( (abs.(ctmp)).^2 ) * pw.Ω / Npoints
        @printf("ibeta, integ_prj = %3d %18.10f %18.10f %18.10f %18.10f %18.10f\n",
                ibeta, integ_prj, integ_prj_im, integ_prj_abs, norm_G.re, norm_G.im)
    end


    Nstates = 4
    Focc = 2.0*ones(Nstates)
    srand(1234)
    psi = rand(Ngwx,Nstates) + im*rand(Ngwx,Nstates)
    psi = ortho_gram_schmidt(psi)

    # calculate < betaNL | psi >
    betaNL_psi = zeros( Complex128, Natoms, Nstates, NbetaNL )
    for ia = 1:Natoms
        for ist = 1:Nstates
            for ibeta = 1:NbetaNL
                betaNL_psi[ia,ist,ibeta] = dot( betaNL[:,ibeta], psi[:,ist] )
            end
        end
    end


    # calculate op_V_NL | psi >
    Vpsi = zeros( Complex128, Ngwx, Nstates )

    for ist = 1:Nstates
        for ia = 1:Natoms
            isp = atm2species[ia]
            for l = 0:psp.lmax
            for m = -l:l
                for iprj = 1:psp.Nproj_l[l+1]
                for jprj = 1:psp.Nproj_l[l+1]
                    ibeta = prj2beta[iprj,ia,l+1,m+psp.lmax+1]
                    jbeta = prj2beta[jprj,ia,l+1,m+psp.lmax+1]
                    hij = psp.h[l+1,iprj,jprj]
                    Vpsi[:,ist] = Vpsi[:,ist] + hij * betaNL[:,ibeta] * betaNL_psi[ia,ist,jbeta]
                end # iprj
                end # jprj
            end # m
            end # l
        end
    end
    println("sum Vpsi = ", sum(Vpsi))

    # calculate E_NL
    E_ps_NL = 0.0
    for ist = 1:Nstates
        enl1 = 0.0
        for ia = 1:Natoms
            isp = atm2species[ia]
            for l = 0:psp.lmax
            for m = -l:l
                for iprj = 1:psp.Nproj_l[l+1]
                for jprj = 1:psp.Nproj_l[l+1]
                    ibeta = prj2beta[iprj,ia,l+1,m+psp.lmax+1]
                    jbeta = prj2beta[jprj,ia,l+1,m+psp.lmax+1]
                    hij = psp.h[l+1,iprj,jprj]
                    enl1 = enl1 + hij*real(conj(betaNL_psi[ia,ist,ibeta])*betaNL_psi[ia,ist,jbeta])
                end
                end
            end # m
            end # l
        end
        E_ps_NL = E_ps_NL + Focc[ist]*enl1
    end
    println("E ps NL = ", E_ps_NL)

end

test_main()
