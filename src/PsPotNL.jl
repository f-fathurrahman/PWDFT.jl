struct PsPotNL
    NbetaNL::Int64
    prj2beta::Array{Int64,4}
    betaNL::Array{ComplexF64,3}
end
# XXX: betaNL should have similar structure with psiks
# Array{Array{ComplexF64,2},1}
# calc_betaNL_psi should utilize matrix multiplication


function PsPotNL()
    # return dummy PsPotNL
    return PsPotNL(0, zeros(Int64,1,1,1,1), zeros(ComplexF64,1,1,1) )
end

function PsPotNL( atoms::Atoms, pw::PWGrid, Pspots::Array{PsPot_GTH,1}; check_norm=false )

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    atpos = atoms.positions

    kpoints = pw.gvecw.kpoints

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
        psp = Pspots[isp]
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
        return PsPotNL(0, zeros(Int64,1,1,1,1), zeros(ComplexF64,1,1,1) )
    end

    Nkpt = kpoints.Nkpt
    k = kpoints.k
    Ngw = pw.gvecw.Ngw
    Ngwx = pw.gvecw.Ngwx
    betaNL = zeros( ComplexF64, Ngwx, NbetaNL, Nkpt )
    G = pw.gvec.G
    g = zeros(3)

    for ik = 1:Nkpt
        ibeta = 0
        for ia = 1:Natoms
            isp = atm2species[ia]
            psp = Pspots[isp]
            for l = 0:psp.lmax
            for iprj = 1:psp.Nproj_l[l+1]
            for m = -l:l
                ibeta = ibeta + 1
                idx_gw2g = pw.gvecw.idx_gw2g[ik]
                for igk = 1:Ngw[ik]
                    ig = idx_gw2g[igk]
                    g[1] = G[1,ig] + k[1,ik]
                    g[2] = G[2,ig] + k[2,ik]
                    g[3] = G[3,ig] + k[3,ik]
                    Gm = norm(g)
                    GX = atpos[1,ia]*g[1] + atpos[2,ia]*g[2] + atpos[3,ia]*g[3]
                    Sf = cos(GX) - im*sin(GX)
                    betaNL[igk,ibeta,ik] =
                    (-1.0*im)^l * Ylm_real(l,m,g)*eval_proj_G(psp,l,iprj,Gm,pw.CellVolume)*Sf
                end
            end
            end
            end
        end
    end  # kpoints

    if check_norm
        check_betaNL_norm( pw, betaNL, kpoints )
    end

    return PsPotNL( NbetaNL, prj2beta, betaNL )

end


function check_betaNL_norm(
    pw::PWGrid,
    betaNL::Array{ComplexF64,3},
    kpoints::KPoints
)

    Npoints = prod(pw.Ns)
    NbetaNL = size(betaNL)[2]

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
            norm_G = dot( betaNL[1:Ngw[ik],ibeta,ik], betaNL[1:Ngw[ik],ibeta,ik] )
            ctmp[:] .= 0.0 + im*0.0
            ctmp[idx_gw2r] = betaNL[1:Ngw[ik],ibeta,ik]
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


function calc_betaNL_psi( ik::Int64, betaNL::Array{ComplexF64,3}, psi::Array{ComplexF64,2} )
    Nstates = size(psi)[2]
    NbetaNL = size(betaNL)[2]
    Ngw_ik = size(psi)[1]
    betaNL_psi = zeros( ComplexF64, Nstates, NbetaNL )
    for ist = 1:Nstates
        for ibeta = 1:NbetaNL
            betaNL_psi[ist,ibeta] = dot( betaNL[1:Ngw_ik,ibeta,ik], psi[:,ist] )
        end
    end
    return betaNL_psi
end

function calc_betaNL_psi( ik::Int64, betaNL::Array{ComplexF64,3}, psi::Array{ComplexF64,1} )
    NbetaNL = size(betaNL)[2]
    Ngw_ik = size(psi)[1]
    betaNL_psi = zeros( ComplexF64, NbetaNL )
    for ibeta = 1:NbetaNL
        betaNL_psi[ibeta] = dot( betaNL[1:Ngw_ik,ibeta,ik], psi )
    end
    return betaNL_psi
end

include("PsPotNL_io.jl")
