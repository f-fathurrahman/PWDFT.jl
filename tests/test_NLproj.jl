using PWDFT

function test_main()
    pw = PWGrid(30.0, 20.0*diagm(ones(3)))
    println(pw)

    # Atoms
    atoms = init_atoms_xyz("Ge.xyz")
    println(atoms)

    # Structure factor
    strf = calc_strfact( atoms, pw )

    zvals = [4.0]
    psp = PsPot_GTH("../pseudopotentials/pade_gth/Ge-q4.gth")
    println(psp)

    Ngwx = pw.gvecw.Ngwx
    idx = pw.gvecw.idx_gw2r
    gwave = pw.gvec.G[:,idx]

    println(size(gwave))

    Nstates = 4
    psi = rand(Ngwx,Nstates) + im*rand(Ngwx,Nstates)
    psi = ortho_gram_schmidt(psi)

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
                        betaNL[ig,ibeta] = Ylm_complex(l,m,g) * eval_proj_G( psp, l, iprj, Gm, pw.Ω ) * Sf
"""                        ylm = Ylm_real(l,m,g)
                        prj = eval_proj_G( psp, l, iprj, Gm, pw.Ω )
                        betaNL[ig,ibeta] = ylm*prj*Sf
                        if ibeta == 4
                            @printf("\n%18.10f %18.10f %18.10f %18.10f\n", g[1], g[2], g[3], Gm)
                            @printf("%18.10f %18.10f %18.10f %18.10f\n", ylm, prj, betaNL[ig,ibeta].re, betaNL[ig,ibeta].im) 
                        end"""
                        #betaNL[ig,ibeta] = Ylm_complex(l,m,g) * exp( -0.5*Gm^2 ) * Sf
                    end
                end
            end
        end
    end

    Npoints = prod(pw.Ns)

    for ibeta = 1:NbetaNL
        ctmp = zeros( Complex128, Npoints )
        ctmp[idx] = betaNL[:,ibeta]
        ctmp = G_to_R(pw, ctmp)*Npoints
        data3d = real(ctmp)
        data3d_im = imag(ctmp)
        
        if ibeta == 5
            LatVecs = pw.LatVecs
            write_xsf( "ATOMS.xsf", LatVecs/ANG2BOHR, atpos/ANG2BOHR, molecule=false, atsymbs=atoms.atsymbs )
            write_xsf_data3d_crystal( "ATOMS.xsf", pw.Ns, LatVecs/ANG2BOHR, data3d )
        end
        #sumBeta = sum(betaNL[:,ibeta])
        #println("ibeta = ", ibeta, " sumBeta = ", sumBeta)

        integ_prj = sum( data3d.^2 ) * pw.Ω / Npoints
        integ_prj_im = sum( data3d_im.^2 ) * pw.Ω / Npoints
        integ_prj_abs = sum( (abs.(ctmp)).^2 ) * pw.Ω / Npoints
        @printf("ibeta, integ_prj = %3d %18.10f %18.10f %18.10f\n",
                ibeta, integ_prj, integ_prj_im, integ_prj_abs)
    end


"""
    l = 2
    m = 0
    iprj = 1
    for ig = 1:Ngwx
        g = gwave[:,ig]
        Gm = norm(g)
        betaNL[ig] = Ylm_real( l, m, g ) * eval_proj_G( psp, l, iprj, Gm, Ω )
        #@printf("%18.10f %18.10f + %18.10fim\n", Gm, betaNL[ig].re, betaNL[ig].im)
    end
"""


    println("rrl = ", psp.rc)
    println("sum betaNL = ", sum(betaNL))


end

test_main()
