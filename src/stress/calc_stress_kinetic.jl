function calc_stress_kinetic!(pw, electrons, psiks, stress_kin)

    G = pw.gvec.G
    idx_gw2g = pw.gvecw.idx_gw2g
    Ngw = pw.gvecw.Ngw
    Ngwx = maximum(Ngw)
    kfac = ones(Float64, Ngwx) # assume qcutz=0
    Gk = zeros(Float64, Ngwx, 3)
    #
    k = pw.gvecw.kpoints.k
    Nkpt = pw.gvecw.kpoints.Nkpt
    wk = pw.gvecw.kpoints.wk
    #
    Nstates = electrons.Nstates
    Nspin = electrons.Nspin
    Focc = electrons.Focc

    fill!(stress_kin, 0.0)

    for ispin in 1:Nspin, ik in 1:Nkpt
        ikspin = ik + (ispin - 1)*Nkpt
        #
        for igw in 1:Ngw[ik]
            ig = idx_gw2g[ik][igw]
            for i in 1:3
                Gk[igw,i] = G[i,ig] + k[i,ik]
            end
            # NOTE: reversed the dimension of Gk[1:3,1:Ngwx] to Gk[1:Ngwx,1:3]
        end
        #@printf("ik, sum(gk) = %4d %18.10f\n", ik, sum(Gk[1:Ngw[ik],:]))
        #
        psi = psiks[ikspin]
        # kinetic contribution
        for l in 1:3, m in 1:l, ist in 1:Nstates, igw in 1:Ngw[ik]
            GlGm = Gk[igw,l] * Gk[igw,m] * kfac[igw]
            psi2 = real( conj(psi[igw,ist]) * psi[igw,ist] )
            stress_kin[l,m] += wk[ik] * GlGm * psi2 * 2 # factor of 2 for spin degeneracy
        end
        #println("wk ik = ", wk[ik]*2)
    end

    for l in 1:3, m in 1:(l-1)
        stress_kin[m,l] = stress_kin[l,m]
    end
    # Scale by 1/pw.CellVolume
    stress_kin[:,:] .*= (1/pw.CellVolume)

    return
end