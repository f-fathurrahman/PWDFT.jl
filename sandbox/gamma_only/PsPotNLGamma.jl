struct PsPotNLGamma
    NbetaNL::Int64
    prj2beta::Array{Int64,4}
    betaNL::Array{ComplexF64,2}
end

function PsPotNLGamma()
    # return dummy PsPotNLGamma
    betaNL = zeros(ComplexF64,1,1)
    return PsPotNLGamma(0, zeros(Int64,1,1,1,1), betaNL )
end

function PsPotNLGamma( atoms::Atoms, pw::PWGridGamma, pspots::Array{PsPot_GTH,1} )

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    atpos = atoms.positions

    # 4: indexed from 0:3
    # 0, 1, 2, 3  -> l indexed
    # 1, 2, 3, 4  -> l + 1

    # -3:3
    # -3, -2, -1, 0, 1, 2, 3  -> m
    #  1,  2,  3, 4, 5, 6, 7  -> 4 + m, lmax = 3 + 1

    # -2, -1, 0, 1, 2  -> m
    #  1,  2, 3, 4, 5  -> 3 + m, lmax = 2 + 1

    prj2beta = Array{Int64}(undef,3,Natoms,4,7)
    prj2beta[:,:,:,:] .= -1   # set to invalid index

    NbetaNL = 0
    for ia = 1:Natoms
        isp = atm2species[ia]
        psp = pspots[isp]
        for l = 0:psp.lmax
            for m = -l:l
                for iprj = 1:psp.Nproj_l[l+1]
                    NbetaNL = NbetaNL + 1
                    prj2beta[iprj,ia,l+1,m+psp.lmax+1] = NbetaNL
                end
            end
        end
    end

    # No nonlocal components
    if NbetaNL == 0
        # return dummy PsPotNL
        betaNL = zeros(ComplexF64,1,1)
        return PsPotNLGamma(0, zeros(Int64,1,1,1,1), betaNL )
    end

    Ngw = pw.gvecw.Ngw
    betaNL = zeros(ComplexF64, Ngw, NbetaNL)

    G = pw.gvec.G
    G2 = pw.gvec.G2
    g = zeros(Float64,3)

    idx_gw2g = pw.gvecw.idx_gw2g
    ibeta = 0

    for ia = 1:Natoms
        isp = atm2species[ia]
        psp = pspots[isp]
        for l = 0:psp.lmax
        for m = -l:l
        for iprj = 1:psp.Nproj_l[l+1]
            ibeta = ibeta + 1
            for igk = 1:Ngw
                ig = idx_gw2g[igk]
                Gm = sqrt(G2[ig])
                # to avoid making slices
                g[1] = G[1,ig]
                g[2] = G[2,ig]
                g[3] = G[3,ig]
                GX = atpos[1,ia]*g[1] + atpos[2,ia]*g[2] + atpos[3,ia]*g[3]
                Sf = cos(GX) - im*sin(GX)
                betaNL[igk,ibeta] =
                (-1.0*im)^l * Ylm_real(l,m,g)*eval_proj_G(psp,l,iprj,Gm,pw.CellVolume)*Sf
            end
        end
        end
        end
    end

    return PsPotNLGamma( NbetaNL, prj2beta, betaNL )

end


function calc_betaNL_psi( pspotNL::PsPotNLGamma, psi::Array{ComplexF64,2} )
    c = psi' * pspotNL.betaNL
    
    Nstates = size(psi,2)
    v1 = zeros(ComplexF64,Nstates)
    for ist in 1:Nstates
        v1[ist] = psi[1,ist]
    end
    
    NbetaNL = size(pspotNL.betaNL,2)
    v2 = zeros(ComplexF64,NbetaNL)
    for i in 1:NbetaNL
        v2[i] = pspotNL.betaNL[1,i]
    end

    return c + conj(c) - v1*v2'
end

