struct PsPotNL
    NbetaNL::Int64
    prj2beta::Array{Int64,4}
    betaNL::Vector{Matrix{ComplexF64}}
    are_ultrasoft::Vector{Bool}
    are_paw::Vector{Bool}
end

function PsPotNL()
    # return dummy PsPotNL
    betaNL = Vector{Matrix{ComplexF64}}(undef,1)
    betaNL[1] = zeros(ComplexF64,1,1)
    return PsPotNL(0, zeros(Int64,1,1,1,1), betaNL, [false], [false] )
end

function PsPotNL( atoms::Atoms, pw::PWGrid, pspots::Vector{PsPot_GTH}; check_norm=false )

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    atpos = atoms.positions
    kpoints = pw.gvecw.kpoints
    are_ultrasoft = zeros(Bool, atoms.Nspecies)
    are_paw = zeros(Bool, atoms.Nspecies)

    NbetaNL, prj2beta = _init_prj2beta(Natoms, atm2species, pspots)

    # No nonlocal components
    if NbetaNL == 0
        # return dummy PsPotNL
        betaNL = Array{Array{ComplexF64,2},1}(undef,1)
        betaNL[1] = zeros(ComplexF64,1,1)
        return PsPotNL(0, zeros(Int64,1,1,1,1), betaNL, are_ultrasoft, are_paw )
    end

    Nkpt = kpoints.Nkpt
    k = kpoints.k
    Ngw = pw.gvecw.Ngw
    Ngwx = pw.gvecw.Ngwx
    betaNL = Array{Array{ComplexF64,2},1}(undef,Nkpt)
    for ik = 1:Nkpt
        betaNL[ik] = zeros(ComplexF64, Ngw[ik], NbetaNL)
    end
    G = pw.gvec.G
    g = zeros(3)

    for ik in 1:Nkpt
        ibeta = 0
        for ia in 1:Natoms
            isp = atm2species[ia]
            psp = pspots[isp]
            for l in 0:psp.lmax, m in -l:l
                for iprj in 1:psp.Nproj_l[l+1]
                    ibeta = ibeta + 1
                    idx_gw2g = pw.gvecw.idx_gw2g[ik]
                    for igk in 1:Ngw[ik]
                        ig = idx_gw2g[igk]
                        g[1] = G[1,ig] + k[1,ik]
                        g[2] = G[2,ig] + k[2,ik]
                        g[3] = G[3,ig] + k[3,ik]
                        Gm = norm(g)
                        GX = atpos[1,ia]*g[1] + atpos[2,ia]*g[2] + atpos[3,ia]*g[3]
                        Sf = cos(GX) - im*sin(GX)
                        betaNL[ik][igk,ibeta] =
                        (-1.0*im)^l * Ylm_real(l,m,g)*eval_proj_G(psp,l,iprj,Gm,pw.CellVolume)*Sf
                    end
                end
            end
        end
    end  # kpoints

    if check_norm
        check_betaNL_norm( pw, betaNL, kpoints )
    end

    return PsPotNL( NbetaNL, prj2beta, betaNL, are_ultrasoft, are_paw )

end


#=
# For PsPot_UPF
# XXX Reduce redundant codes
function PsPotNL( atoms::Atoms, pw::PWGrid, pspots::Array{PsPot_UPF,1}; check_norm=false )

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    atpos = atoms.positions
    kpoints = pw.gvecw.kpoints

    NbetaNL, prj2beta = _init_prj2beta(Natoms, atm2species, pspots)

    # No nonlocal components
    if NbetaNL == 0
        # return dummy PsPotNL
        betaNL = Array{Array{ComplexF64,2},1}(undef,1)
        betaNL[1] = zeros(ComplexF64,1,1)
        return PsPotNL(0, zeros(Int64,1,1,1,1), betaNL )
    end

    Nkpt = kpoints.Nkpt
    k = kpoints.k
    Ngw = pw.gvecw.Ngw
    Ngwx = pw.gvecw.Ngwx
    betaNL = Array{Array{ComplexF64,2},1}(undef,Nkpt)
    for ik = 1:Nkpt
        betaNL[ik] = zeros(ComplexF64, Ngw[ik], NbetaNL)
    end
    G = pw.gvec.G
    g = zeros(3)

    for ik in 1:Nkpt
        idx_gw2g = pw.gvecw.idx_gw2g[ik]
        for ia in 1:Natoms
            isp = atm2species[ia]
            psp = pspots[isp]
            iprjl = 0
            # FIXME: This loop is not optimal, probably should loop over combined index lm
            for l in 0:psp.lmax, iprj in 1:psp.Nproj_l[l+1]
                iprjl = iprjl + 1 # increment projector index
                for m in -l:l
                    ibeta = prj2beta[iprj,ia,l+1,m+psp.lmax+1]
                    for igk in 1:Ngw[ik]
                        ig = idx_gw2g[igk]
                        g[1] = G[1,ig] + k[1,ik]
                        g[2] = G[2,ig] + k[2,ik]
                        g[3] = G[3,ig] + k[3,ik]
                        Gm = norm(g)
                        GX = atpos[1,ia]*g[1] + atpos[2,ia]*g[2] + atpos[3,ia]*g[3]
                        Sf = cos(GX) - im*sin(GX)
                        betaNL[ik][igk,ibeta] =
                        (-1.0*im)^l * Ylm_real(l,m,g)*eval_proj_G(psp, iprjl, Gm)*Sf
                    end
                end
            end
        end
    end

    if check_norm
        check_betaNL_norm( pw, betaNL, kpoints )
    end

    return PsPotNL( NbetaNL, prj2beta, betaNL )

end
=#

function _init_prj2beta(Natoms::Int64, atm2species, pspots)
    # 4: indexed from 0:3
    # 0, 1, 2, 3  -> l indexed
    # 1, 2, 3, 4  -> l + 1

    # -3:3
    # -3, -2, -1, 0, 1, 2, 3  -> m
    #  1,  2,  3, 4, 5, 6, 7  -> 4 + m, lmax = 3 + 1

    # -2, -1, 0, 1, 2  -> m
    #  1,  2, 3, 4, 5  -> 3 + m, lmax = 2 + 1

    # Fix for full-relativistic ONCV
    # FIXME: Note that full-relativistic pspot is not yet supported
    nprojlmax = 1
    for psp in pspots
        nn = maximum(psp.Nproj_l)
        if nprojlmax < nn
            nprojlmax = nn
        end
    end
    prj2beta = Array{Int64}(undef,nprojlmax,Natoms,4,7)
    fill!(prj2beta, -1)   # set to invalid index

    NbetaNL = 0
    for ia in 1:Natoms
        isp = atm2species[ia]
        psp = pspots[isp]
        for l in 0:psp.lmax, m in -l:l
            for iprj in 1:psp.Nproj_l[l+1]
                NbetaNL = NbetaNL + 1
                prj2beta[iprj,ia,l+1,m+psp.lmax+1] = NbetaNL
            end
        end
    end

    #for ia in 1:Natoms
    #    isp = atm2species[ia]
    #    psp = pspots[isp]
    #    for l in 0:psp.lmax, iprj in 1:psp.Nproj_l[l+1]
    #        for m in -l:l
    #            NbetaNL = NbetaNL + 1
    #            prj2beta[iprj,ia,l+1,m+psp.lmax+1] = NbetaNL
    #        end
    #    end
    #end

    return NbetaNL, prj2beta
end



function check_betaNL_norm(
    pw::PWGrid,
    betaNL::Array{Array{ComplexF64,2},1},
    kpoints::KPoints
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


function calc_betaNL_psi(
    ik::Int64,
    betaNL::Array{Array{ComplexF64,2},1},
    psi::AbstractArray{ComplexF64}
)

    Nstates = size(psi,2)
    NbetaNL = size(betaNL[1],2)

    betaNL_psi = zeros( ComplexF64, Nstates, NbetaNL )
    betaNL_psi[:,:] = conj( psi' * betaNL[ik] ) # we calculate < betaNL | psi |
    return betaNL_psi
end

import Base: show
function show( io::IO, pspotNL::PsPotNL )
    
    println("--------")
    println("PsPotNL:")
    println("--------")
    
    println("NbetaNL  = ", pspotNL.NbetaNL)

    return
end