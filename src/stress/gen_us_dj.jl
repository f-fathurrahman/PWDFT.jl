function _gen_us_dj!(ik, atoms, pw, pspots, pspotNL, dvkb)

    dq = 0.01 # HARDCODED

    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atpos = atoms.positions
    atm2species = atoms.atm2species

    G = pw.gvec.G
    idx_gw2g = pw.gvecw.idx_gw2g
    Ngwk = pw.gvecw.Ngw[ik]
    k = pw.gvecw.kpoints.k

    nh = pspotNL.nh
    indv = pspotNL.indv
    nhtol = pspotNL.nhtol
    nhtolm = pspotNL.nhtolm
    lmaxkb = pspotNL.lmaxkb
    prj_interp_tables = pspotNL.prj_interp_tables

    Gk = zeros(Float64, 3, Ngwk)    
    Gk2 = zeros(Float64, Ngwk)
    for igk in 1:Ngwk
        ig = idx_gw2g[ik][igk]
        Gk[1,igk] = k[1,ik] + G[1,ig]
        Gk[2,igk] = k[2,ik] + G[2,ig]
        Gk[3,igk] = k[3,ik] + G[3,ig]
        Gk2[igk] = Gk[1,igk]^2 + Gk[2,igk]^2 + Gk[3,igk]^2
    end

    ylm = zeros(Float64, Ngwk, (lmaxkb+1)^2)
    Ylm_real_qe!(lmaxkb, Gk, ylm) # Ylm_real_qe accept l value starting from 0


    djl = Vector{Matrix{Float64}}(undef, Nspecies)
    for isp in 1:Nspecies
        Nproj = pspots[isp].Nproj
        djl[isp] = zeros(Float64, Ngwk, Nproj)
    end
    #
    for isp in 1:Nspecies
        psp = pspots[isp]
        tab = prj_interp_tables[isp]
        for ibeta in 1:pspots[isp].Nproj, igk in 1:Ngwk
            Gm = sqrt(Gk2[igk])
            px = Gm/dq - floor(Int64, Gm/dq)
            ux = 1.0 - px
            vx = 2.0 - px
            wx = 3.0 - px
            i0 = floor(Int64, Gm/dq) + 1
            i1 = i0 + 1
            i2 = i0 + 2
            i3 = i0 + 3
            djl[isp][igk,ibeta] = ( tab[i0,ibeta] * ( -vx*wx - ux*wx - ux*vx )/6 +
                                    tab[i1,ibeta] * (  vx*wx - px*wx - px*vx )/2 -
                                    tab[i2,ibeta] * (  ux*wx - px*wx - px*ux )/2 +
                                    tab[i3,ibeta] * (  ux*vx - px*vx - px*ux )/6 )/dq
        end
    end


    Sf = zeros(ComplexF64, Ngwk)
    ikb = 0
    fill!(dvkb, 0.0 + im*0.0)
    for isp in 1:Nspecies, ia in 1:Natoms
        #
        # skip if current atom is not of species isp
        if atm2species[ia] != isp
            continue
        end
        #
        for igk in 1:Ngwk
            GkX = atpos[1,ia]*Gk[1,igk] + atpos[2,ia]*Gk[2,igk] + atpos[3,ia]*Gk[3,igk]
            Sf[igk] = cos(GkX) - im*sin(GkX)
        end
        #
        for ih in 1:nh[isp]
            ikb = ikb + 1
            ibeta = indv[ih,isp]
            l = nhtol[ih,isp]
            lm = nhtolm[ih,isp]
            pref = (-im)^l
            for igk in 1:Ngwk
                dvkb[igk,ikb] = djl[isp][igk,ibeta] * Sf[igk] * ylm[igk,lm] * pref
            end
        end
    end


    return

end