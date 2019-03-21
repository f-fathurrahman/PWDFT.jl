function PsPotNL_v2(
    atoms::Atoms,
    pw::PWGrid,
    pspots::Array{PsPot_GTH,1};
    check_norm=false
)
    kpoints = pw.gvecw.kpoints
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    atpos = atoms.positions

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
                    @printf("ibeta=%3d iprj=%3d l=%2d m=%2d\n", NbetaNL, iprj, l, m)
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
            psp = pspots[isp]
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
                    Ylm_real(l,m,g) * eval_proj_G(psp,l,iprj,Gm,pw.CellVolume)*Sf
                end
            end
            end
            end
        end
    end  # kpoints

    if check_norm
        @time check_betaNL_norm( pw, kpoints, betaNL )
        @time check_betaNL_norm( pw, kpoints, betaNL )
    end

    return PsPotNL( NbetaNL, prj2beta, betaNL )

end

function check_betaNL_norm(
    pw::PWGrid,
    kpoints::KPoints,
    betaNL::Array{ComplexF64,3}
)

    Npoints = prod(pw.Ns)
    NbetaNL = size(betaNL)[2]

    Nkpt = kpoints.Nkpt
    Ngw = pw.gvecw.Ngw

    dVol = pw.CellVolume/Npoints
    println("dVol = ", dVol)
    #
    # Check normalization in real space and reciprocal space
    #
    ctmp = zeros(ComplexF64, Npoints)

    for ik = 1:Nkpt
        #
        idx_gw2r = pw.gvecw.idx_gw2r[ik]
        #
        for ibeta = 1:NbetaNL
            
            f = @view betaNL[:,ibeta,ik]
            norm_G = dot(f,f)
            
            ctmp[:] .= 0.0 + im*0.0
            ctmp[idx_gw2r] = betaNL[:,ibeta,ik]
            
            # This is required for normalization for wavefunction-like quantity
            ctmp[:] = G_to_R(pw, ctmp)*Npoints/sqrt(pw.CellVolume)
                        
            integ_prj = 0.0
            integ_prj_im = 0.0
            integ_prj_abs = 0.0
            for ip = 1:Npoints
                c = ctmp[ip]
                integ_prj = integ_prj + c.re*c.re
                integ_prj_im = integ_prj_im + c.im*c.im
                integ_prj_abs = integ_prj_abs + c.re*c.re + c.im*c.im
            end

            integ_prj = integ_prj*dVol
            integ_prj_im = integ_prj_im*dVol
            integ_prj_abs = integ_prj_abs*dVol
            
            #integ_prj = sum( data3d.^2 ) * pw.CellVolume / Npoints
            #integ_prj_im = sum( data3d_im.^2 ) * pw.CellVolume / Npoints
            #integ_prj_abs = sum( (abs.(ctmp)).^2 ) * pw.CellVolume / Npoints
            
            @printf("%3d %18.10f %18.10f %18.10f %18.10f %18.10f\n",
                     ibeta, integ_prj, integ_prj_im, integ_prj_abs, norm_G.re, norm_G.im)
        end
    end

end