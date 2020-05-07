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

function PsPotNLGamma( atoms::Atoms, pw::PWGridGamma, pspots::Array{PsPot_GTH,1}; check_norm=false )

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
            for iprj = 1:psp.Nproj_l[l+1]
                for m = -l:l
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

    #betaNL_full = zeros(ComplexF64,2*Ngw-1,NbetaNL)

    for ia = 1:Natoms
        isp = atm2species[ia]
        psp = pspots[isp]
        for l = 0:psp.lmax
        for iprj = 1:psp.Nproj_l[l+1]
        for m = -l:l
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

                #betaNL_full[igk,ibeta] =
                #(-1.0*im)^l * Ylm_real(l,m,g)*eval_proj_G(psp,l,iprj,Gm,pw.CellVolume)*Sf
                #if igk != 1
                #    g[1] = -G[1,ig]
                #    g[2] = -G[2,ig]
                #    g[3] = -G[3,ig]
                #    GX = atpos[1,ia]*g[1] + atpos[2,ia]*g[2] + atpos[3,ia]*g[3]
                #    Sf = cos(GX) - im*sin(GX)
                #    betaNL_full[Ngw+igk-1,ibeta] =
                #    (-1.0*im)^l * Ylm_real(l,m,g)*eval_proj_G(psp,l,iprj,Gm,pw.CellVolume)*Sf
                #end
            end
        end
        end
        end
    end

    #if check_norm
    #    check_betaNL_norm( pw, betaNL, kpoints )
    #end

    #println("sum betaNL      = ", sum(betaNL))
    #println("sum betaNL_full = ", sum(betaNL_full))

    return PsPotNLGamma( NbetaNL, prj2beta, betaNL )

end

#=
function check_betaNL_norm(
    pw::PWGrid,
    betaNL::Array{ComplexF64,2}
)

    Npoints = prod(pw.Ns)
    NbetaNL = size(betaNL[1],2)

    Nkpt = kpoints.Nkpt
    Ngw = pw.gvecw.Ngw

    #
    # Check normalization in real space and reciprocal space
    #
    ctmp = zeros( ComplexF64, Npoints )
    for ik = 1:Nkpt
        #
        idx_gw2r = pw.gvecw.idx_gw2r[ik]
        #
        for ibeta = 1:NbetaNL
            
            norm_G = dot( betaNL[ik][:,ibeta], betaNL[ik][:,ibeta] )
            
            ctmp[:] .= 0.0 + im*0.0
            
            ctmp[idx_gw2r] = betaNL[ik][:,ibeta]

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
    end

end
=#

#=
function calc_betaNL_psi(
    betaNL::Array{ComplexF64,2},
    psi::Array{ComplexF64,2}
)

    Nstates = size(psi)[2]
    NbetaNL = size(betaNL[1],2)

    betaNL_psi = zeros( ComplexF64, Nstates, NbetaNL )
    betaNL_psi[:,:] = conj( psi' * betaNL[ik] )
    return betaNL_psi
end

function calc_betaNL_psi(
    ik::Int64,
    betaNL::Array{Array{ComplexF64,2},1},
    psi::Array{ComplexF64,1}
)
    NbetaNL = size(betaNL[1],2)

    betaNL_psi = zeros( ComplexF64, NbetaNL )
    betaNL_psi[:] = conj( psi' * betaNL[ik] )
    return betaNL_psi
end
=#
