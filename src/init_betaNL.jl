function init_betaNL( pw::PWGrid, Pspots::Array{PsPot_GTH,1}, atoms::Atoms, NbetaNL::Int64 )
    
    Ngwx = pw.gvecw.Ngwx
    idx = pw.gvecw.idx_gw2r
    gwave = pw.gvec.G[:,idx]

    atm2species = atoms.atm2species
    atpos = atoms.atpos

    betaNL = zeros( Complex128, Ngwx, NbetaNL )

    ibeta = 0
    for ia = 1:Natoms
        isp = atm2species[ia]
        psp = Pspots[isp]
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

    return betaNL

end
