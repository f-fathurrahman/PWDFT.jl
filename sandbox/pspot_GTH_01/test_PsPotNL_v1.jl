using LinearAlgebra
using Random
using Printf
using PWDFT

#=
Test for handling of PsPotNL.
The new implementation in Hamiltonian is slightly more efficient than this.
=#
function test_main()
    LatVecs = gen_lattice_sc(20.0)
    pw = PWGrid(80.0, LatVecs)
    println(pw)

    # Atoms
    #atoms = init_atoms_xyz("../structures/PtO.xyz")
    atoms = init_atoms_xyz_string("""
        1

        Cs  0.0  0.0  0.0
        """)
    atoms.LatVecs = LatVecs
    println(atoms)

    # Structure factor
    strf = calc_strfact( atoms, pw )

    Nspecies = atoms.Nspecies
    PsPots = Array{PsPot_GTH}(undef,Nspecies)

    # XXX Need to match Nspecies
    #PsPots[1] = PsPot_GTH("../pseudopotentials/pade_gth/Pt-q10.gth")
    #PsPots[2] = PsPot_GTH("../pseudopotentials/pade_gth/O-q6.gth")

    #PsPots[1] = PsPot_GTH("../pseudopotentials/pade_gth/Ba-q10.gth")

    PsPots[1] = PsPot_GTH("../pseudopotentials/pade_gth/Cs-q9.gth")

    for psp in PsPots
        println(psp)
    end

    ik = 1
    Ngwx = pw.gvecw.Ngwx
    idx_gw2g = pw.gvecw.idx_gw2g[ik]
    idx_gw2r = pw.gvecw.idx_gw2r[ik]
    gwave = pw.gvec.G[:,idx_gw2g]

    Natoms = atoms.Natoms

    # 4: indexed from 0:3
    # 0, 1, 2, 3  -> l indexed
    # 1, 2, 3, 4  -> l + 1

    # -3:3
    # -3, -2, -1, 0, 1, 2, 3  -> m
    #  1,  2,  3, 4, 5, 6, 7  -> 4 + m, lmax = 3 + 1

    # -2, -1, 0, 1, 2  -> m
    #  1,  2, 3, 4, 5  -> 3 + m, lmax = 2 + 1

    prj2beta = Array{Int64}( undef, 3, Natoms, 4, 7 )
    prj2beta[:,:,:,:] .= -1   # set to invalid index

    atm2species = atoms.atm2species
    atpos = atoms.positions

    NbetaNL = 0
    for ia = 1:Natoms
        isp = atm2species[ia]
        psp = PsPots[isp]
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

    betaNL = zeros( ComplexF64, Ngwx, NbetaNL )

    ibeta = 0
    for ia = 1:Natoms
        isp = atm2species[ia]
        psp = PsPots[isp]
        for l = 0:psp.lmax
            for iprj = 1:psp.Nproj_l[l+1]
                for m = -l:l
                    ibeta = ibeta + 1
                    for ig = 1:Ngwx
                        g = gwave[:,ig]
                        Gm = norm(g)
                        GX = atpos[1,ia]*g[1] + atpos[2,ia]*g[2] + atpos[3,ia]*g[3]
                        Sf = cos(GX) - im*sin(GX)
                        betaNL[ig,ibeta] = (1.0*im)^l * Ylm_real(l,m,g) * eval_proj_G( psp, l, iprj, Gm, pw.CellVolume ) * Sf
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
        ctmp = zeros( ComplexF64, Npoints )
        ctmp[idx_gw2r] = betaNL[:,ibeta]
        ctmp = G_to_R(pw, ctmp)*Npoints
        data3d = real(ctmp)
        data3d_im = imag(ctmp)
        #
        integ_prj = sum( data3d.^2 ) * pw.CellVolume / Npoints
        integ_prj_im = sum( data3d_im.^2 ) * pw.CellVolume / Npoints
        integ_prj_abs = sum( (abs.(ctmp)).^2 ) * pw.CellVolume / Npoints
        @printf("ibeta, integ_prj = %3d %18.10f %18.10f %18.10f %18.10f %18.10f\n",
                ibeta, integ_prj, integ_prj_im, integ_prj_abs, norm_G.re, norm_G.im)
    end


    Nstates = 4
    Focc = 2.0*ones(Nstates)
    Random.seed!(1234)
    psi = rand(ComplexF64,Ngwx,Nstates)
    psi = ortho_gram_schmidt(psi)

    # calculate < betaNL | psi >
    betaNL_psi = zeros( ComplexF64, Natoms, Nstates, NbetaNL )
    for ia = 1:Natoms
        for ist = 1:Nstates
            for ibeta = 1:NbetaNL
                betaNL_psi[ia,ist,ibeta] = dot( betaNL[:,ibeta], psi[:,ist] )
            end
        end
    end


    # calculate op_V_NL | psi >
    Vpsi = zeros( ComplexF64, Ngwx, Nstates )

    for ist = 1:Nstates
        for ia = 1:Natoms
            isp = atm2species[ia]
            psp = PsPots[isp]
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


    # Build array hij
    hij = zeros(Natoms,NbetaNL,NbetaNL)
    for ia = 1:Natoms
        isp = atm2species[ia]
        psp = PsPots[isp]
        for l = 0:psp.lmax
        for m = -l:l
            for iprj = 1:psp.Nproj_l[l+1]
            for jprj = 1:psp.Nproj_l[l+1]
                ibeta = prj2beta[iprj,ia,l+1,m+psp.lmax+1]
                jbeta = prj2beta[jprj,ia,l+1,m+psp.lmax+1]
                hij[ia,ibeta,jbeta] = psp.h[l+1,iprj,jprj]
            end
            end
        end # m
        end # l
    end

    # calculate E_NL
    E_ps_NL = 0.0
    for ist = 1:Nstates
        enl1 = 0.0
        for ia = 1:Natoms
            isp = atm2species[ia]
            psp = PsPots[isp]
            for l = 0:psp.lmax
            for m = -l:l
                for iprj = 1:psp.Nproj_l[l+1]
                for jprj = 1:psp.Nproj_l[l+1]
                    ibeta = prj2beta[iprj,ia,l+1,m+psp.lmax+1]
                    jbeta = prj2beta[jprj,ia,l+1,m+psp.lmax+1]
                    #hij = psp.h[l+1,iprj,jprj]
                    enl1 = enl1 + hij[ia,ibeta,jbeta]*real(conj(betaNL_psi[ia,ist,ibeta])*betaNL_psi[ia,ist,jbeta])
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
