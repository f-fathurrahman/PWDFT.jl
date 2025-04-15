# Experimenting with USPP

mutable struct PsPotNL_UPF
    lmaxx::Int64
    lqmax::Int64
    lmaxkb::Int64
    nh::Vector{Int64}
    nhm::Int64
    NbetaNL::Int64  # it is originally nkb
    prj_interp_tables::Vector{Matrix{Float64}}
    ap::Array{Float64,3}
    lpx::Array{Int64,2}
    lpl::Array{Int64,3}
    indv::Array{Int64,2}
    nhtol::Array{Int64,2}
    nhtolm::Array{Int64,2}
    ijtoh::Array{Int64,3}
    indv_ijkb0::Vector{Int64}
    Dvan::Array{Float64,3}
    Deeq::Array{Float64,4}
    qradG::Union{Vector{Array{Float64,3}},Nothing}
    qq_nt::Union{Array{Float64,3},Nothing}
    qq_at::Union{Array{Float64,3},Nothing}
    betaNL::Vector{Matrix{ComplexF64}}
    are_ultrasoft::Vector{Bool} # XXX rename to are_uspp?
    are_paw::Vector{Bool}
    becsum::Array{Float64,3}
    paw::Union{Nothing,PAWVariables}
end


function PsPotNL_UPF(
    atoms::Atoms,
    pw::PWGrid,
    pspots::Vector{PsPot_UPF};
    is_gga=false,
    Nspin=1
)

    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species

    # Those are parameters (HARDCODED)
    lmaxx  = 3         # max non local angular momentum (l=0 to lmaxx)
    lqmax = 2*lmaxx + 1  # max number of angular momenta of Q

    # calculate the number of beta functions for each atomic type
    nh = zeros(Int64, Nspecies)
    lmaxkb = -1
    for isp in 1:Nspecies
        psp = pspots[isp]
        for nb in 1:psp.Nproj
            nh[isp] = nh[isp] + 2 * psp.proj_l[nb] + 1
            lmaxkb = max(lmaxkb, psp.proj_l[nb])
        end
    end

    NbetaNL = 0
    for ia in 1:Natoms
       isp = atm2species[ia]
       NbetaNL = NbetaNL + nh[isp]
    end

    # Check for UPF GTH_num pseudopotentials
    #@assert NbetaNL > 0
    # Can still run if NbetaNL == 0
    # Various arrays may be zeros due to NbetaNL == 0

    # Nonlocal projector interpolation tables for each species
    prj_interp_tables = Vector{Matrix{Float64}}(undef,Nspecies)
    for isp in 1:Nspecies
        prj_interp_tables[isp] = _build_prj_interp_table(pspots[isp], pw)
    end

    ap, lpx, lpl = _calc_clebsch_gordan(lmaxkb)

    # calculate the maximum number of beta functions
    nhm = maximum(nh)
    # Some helper indices, for mapping various stuffs
    # Some of these can be made into a jagged array (vector of vector)
    indv = zeros(Int64, nhm, Nspecies)
    nhtol = zeros(Int64, nhm, Nspecies)
    nhtolm = zeros(Int64, nhm, Nspecies)
    indv_ijkb0 = zeros(Int64, Natoms)
    ijtoh = zeros(Int64, nhm, nhm, Nspecies)
    Dvan = zeros(Float64, nhm, nhm, Nspecies)

    ijkb0 = 0
    for isp in 1:Nspecies
        ih = 1
        psp = pspots[isp]
        for nb in 1:psp.Nproj
            l = psp.proj_l[nb]
            for m in 1:(2*l+1)
                nhtol[ih,isp] = l
                nhtolm[ih,isp] = l*l + m
                indv[ih,isp] = nb
                ih = ih + 1
            end
        end
        #
        # ijtoh map augmentation channel indexes ih and jh to composite
        # "triangular" index ijh
        @views ijtoh[:,:,isp] .= -1
        ijv = 0
        for ih in 1:nh[isp], jh in ih:nh[isp]
             ijv = ijv + 1
             ijtoh[ih,jh,isp] = ijv
             ijtoh[jh,ih,isp] = ijv
        end
        #
        # ijkb0 points to the last beta "in the solid" for atom ia-1
        # i.e. ijkb0+1,.. ijkb0+nh(ityp(ia)) are the nh beta functions of
        #      atom ia in the global list of beta functions (ijkb0=0 for ia=1)
        for ia in 1:Natoms
            if atm2species[ia] != isp
                continue
            end
            indv_ijkb0[ia] = ijkb0
            ijkb0 = ijkb0 + nh[isp]
        end

        for ih in 1:nh[isp], jh in 1:nh[isp]
            cond_l  = nhtol[ih,isp] == nhtol[jh,isp]
            cond_lm = nhtolm[ih,isp] == nhtolm[jh,isp]
            if cond_l && cond_lm
                ihs = indv[ih,isp]
                jhs = indv[jh,isp]
                Dvan[ih,jh,isp] = psp.Dion[ihs,jhs]
            end
        end
    end
    # TODO: Extract lm -> (l,m)

    are_ultrasoft = zeros(Bool,Nspecies)
    for isp in 1:Nspecies
        are_ultrasoft[isp] = pspots[isp].is_ultrasoft
    end
    # XXX: should be any_uspp ?
    if all(are_ultrasoft)
        qradG, qq_nt, qq_at = _prepare_aug_charges(
            atoms, pw, pspots, lmaxkb, nhm, nh, indv, nhtolm, lpl, lpx, ap
        )
    else
        # XXX Use zeros or nothing?
        # XXX If we use zeros we need information about the sizes of these arrays
        qradG = nothing
        qq_nt = nothing
        qq_at = nothing
    end

    Nkpt = pw.gvecw.kpoints.Nkpt
    betaNL = Vector{Matrix{ComplexF64}}(undef,Nkpt)
    for ik in 1:Nkpt
        betaNL[ik] = zeros(ComplexF64, pw.gvecw.Ngw[ik], NbetaNL)
        _init_Vnl_KB!(
            ik, atoms, pw, pspots, prj_interp_tables,
            lmaxkb, nh, nhm, nhtol, nhtolm, indv,
            betaNL[ik]
        )
    end

    # The D-matrix
    # Contains NL pspot coefficients with contribution from
    # integral of augmentation charges times Veff)
    # Depends on spin
    Deeq = zeros(Float64, nhm, nhm, Natoms, Nspin)
    # Set to Dvan if no ultrasoft (need this?)
    #if all(.!are_ultrasoft)
        for ia in 1:Natoms
            isp = atm2species[ia]
            nht = nh[isp]
            for ispin in 1:Nspin
                @views Deeq[1:nht,1:nht,ia,ispin] = Dvan[1:nht,1:nht,isp]
            end
        end
    #end


    are_paw = zeros(Bool,Nspecies)
    for isp in 1:Nspecies
        are_paw[isp] = pspots[isp].is_paw
    end
    # FIXME: We do not consider mixing PAW and USPP
    if all(are_paw)
        paw = PAWVariables(atoms, pspots, nhm, is_gga=is_gga, Nspin=Nspin)
    else
        paw = nothing
    end

    Nbecsum = Int64( nhm * (nhm + 1)/2 )
    becsum = zeros(Float64, Nbecsum, Natoms, Nspin)

    return PsPotNL_UPF(
        lmaxx, lqmax, lmaxkb,
        nh, nhm, NbetaNL, prj_interp_tables,
        ap, lpx, lpl,
        indv, nhtol, nhtolm, ijtoh, indv_ijkb0,
        Dvan, Deeq, qradG, qq_nt, qq_at,
        betaNL,
        are_ultrasoft, are_paw,
        becsum, paw
    )

end


function _build_prj_interp_table( psp::PsPot_UPF, pw::PWGrid )

    ecutwfc = pw.ecutwfc
    CellVolume = pw.CellVolume

    cell_factor = 1.0 # XXX HARDCODED
    dq = 0.01 # XXX HARDCODED

    ndm = psp.kkbeta
    Nproj = psp.Nproj

    nqx = floor( Int64, (sqrt(2*ecutwfc)/dq + 4)*cell_factor )

    prj_interp_table = zeros(Float64, nqx, Nproj)

    aux = zeros(Float64, ndm)
    pref = 4*pi/sqrt(CellVolume)

    for ibeta in 1:Nproj
        l = psp.proj_l[ibeta]
        for iq in 1:nqx
            qi = (iq - 1) * dq
            for ir in 1:psp.kkbeta
                jlqr = sphericalbesselj(l, qi*psp.r[ir])
                aux[ir] = psp.proj_func[ir,ibeta] * psp.r[ir] * jlqr
            end
            vqint = integ_simpson( psp.kkbeta, aux, psp.rab )
            prj_interp_table[iq, ibeta] = vqint * pref
        end
    end

    return prj_interp_table
end


# XXX: Need this?
function _init_Vnl_KB!(
    ik::Int64,
    atoms::Atoms,
    pw::PWGrid,
    pspots::Vector{PsPot_UPF},
    pspotNL::PsPotNL_UPF,
    Vnl_KB::Array{ComplexF64,2}
)
    _init_Vnl_KB!(
        ik, atoms, pw, pspots, pspotNL.prj_interp_tables,
        pspotNL.lmaxkb, pspotNL.nh, pspotNL.nhm,
        pspotNL.nhtol, pspotNL.nhtolm, pspotNL.indv,
        Vnl_KB
    )
    return
end



# From init_us_2
function _init_Vnl_KB!(
    ik::Int64,
    atoms::Atoms,
    pw::PWGrid,
    pspots::Vector{PsPot_UPF},
    prj_interp_tables,
    lmaxkb, nh, nhm, nhtol, nhtolm, indv,
    Vnl_KB::Array{ComplexF64,2}
)

    Ngw = pw.gvecw.Ngw
    @assert size(Vnl_KB,1) == Ngw[ik]

    idx_gw2g = pw.gvecw.idx_gw2g
    k = pw.gvecw.kpoints.k
    G = pw.gvec.G

    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atpos = atoms.positions
    atm2species = atoms.atm2species

    Gk = zeros(Float64, 3, Ngw[ik])
    Gk2 = zeros(Float64, Ngw[ik])
    # 
    # igk: index for Gk
    # ig: index for G
    for igk in 1:Ngw[ik]
        ig = idx_gw2g[ik][igk] # index of Gk in G
        Gk[1,igk] = G[1,ig] + k[1,ik]
        Gk[2,igk] = G[2,ig] + k[2,ik]
        Gk[3,igk] = G[3,ig] + k[3,ik]
        Gk2[igk] = Gk[1,igk]^2 +  Gk[2,igk]^2 + Gk[3,igk]^2
        # need Gk2? it is also calculated in Ylm_real_qe!
    end

    ylm = zeros(Float64, Ngw[ik], (lmaxkb+1)^2)
    # Ylm_real_qe accept l value starting from 0 (the actual 'physics' angular momentum number)
    Ylm_real_qe!(lmaxkb, Gk, ylm)


    vq = zeros(Float64, Ngw[ik])
    vkb1 = zeros(Float64, Ngw[ik], nhm)
    Sf = zeros(ComplexF64, Ngw[ik])
    #
    # |beta_lm(q)> = (4pi/omega).Y_lm(q).f_l(q).(i^l).S(q)
    jkb = 0 # index of beta functions
    dq = 0.01 # HARDCODED
    #    
    for isp in 1:Nspecies
        #
        psp = pspots[isp]
        tab = prj_interp_tables[isp]
        #
        # calculate beta in G-space using an interpolation table:
        # f_l(q)=\int _0 ^\infty dr r^2 f_l(r) j_l(q.r)
        for ibeta in 1:psp.Nproj
            #
            for igk in 1:Ngw[ik]
                Gm = sqrt(Gk2[igk])
                px = Gm/dq - floor(Int64, Gm/dq )
                ux = 1.0 - px
                vx = 2.0 - px
                wx = 3.0 - px
                i0 = floor(Int64, Gm/dq ) + 1
                i1 = i0 + 1
                i2 = i0 + 2
                i3 = i0 + 3
                vq[igk] = tab[i0,ibeta] * ux * vx * wx / 6.0 +
                          tab[i1,ibeta] * px * vx * wx / 2.0 -
                          tab[i2,ibeta] * px * ux * wx / 2.0 +
                          tab[i3,ibeta] * px * ux * vx / 6.0
            end

      
            # add spherical harmonic part  (Y_lm(q)*f_l(q)) 
            for ih in 1:nh[isp]
                if ibeta == indv[ih,isp]        
                    lm = nhtolm[ih,isp]
                    for igk in 1:Ngw[ik]
                        vkb1[igk,ih] = ylm[igk,lm] * vq[igk]
                    end
                end
            end

        end

        # vkb1 contains all betas including angular part for type nt
        # now add the structure factor and factor (-i)^l

        # ordering: first all betas for atoms of type 1
        #           then  all betas for atoms of type 2  and so on
        #
        for ia in 1:Natoms
            # skip if this is not the current species index
            if atm2species[ia] != isp
                continue
            end
            #
            for igk in 1:Ngw[ik]
                # XXX use dot product here?
                GkX = atpos[1,ia]*Gk[1,igk] + atpos[2,ia]*Gk[2,igk] + atpos[3,ia]*Gk[3,igk]
                Sf[igk] = cos(GkX) - im*sin(GkX)
            end
            
            for ih in 1:nh[isp]
                # XXXX idx of KB projectors increment here ...
                jkb = jkb + 1
                # No need to offset nhtol by 1. It is already physics's l (start from l=0)
                pref = (-im)^nhtol[ih,isp]
                for igk in 1:Ngw[ik]
                    Vnl_KB[igk,jkb] = vkb1[igk,ih] * Sf[igk] * pref
                end
            end
        end
    
    end

    return

end


# XXX Need to check the unit. It should be the same as chgden.
function _prepare_aug_charges(
    atoms::Atoms, pw::PWGrid, pspots::Vector{PsPot_UPF},
    lmaxkb, nhm, nh, indv, nhtolm, lpl, lpx, ap
)
    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species

    qradG = _calc_qradG(pw, pspots)

    # Compute the qq coefficients by integrating the Q.
    # The qq are the g=0 components of Q

    # FIXME: use Vector{Matrix} instead of 3d array
    qq_nt = zeros(Float64, nhm, nhm, Nspecies) # qq_nt, need qvan2
    qq_at = zeros(Float64, nhm, nhm, Natoms)

    G0 = zeros(3,1) # Needs to be a two dimensional array
    lmaxq = 2*lmaxkb + 1 # using 1-indexing
    ylmk0 = zeros(Float64, 1, lmaxq*lmaxq)
    _lmax = lmaxq - 1 # or 2*lmaxkb
    Ylm_real_qe!(_lmax, G0, ylmk0) # Ylm_real_qe accept l value starting from 0
    qgm = zeros(ComplexF64, 1)
    for isp in 1:Nspecies
        for ih in 1:nh[isp], jh in ih:nh[isp]
            qvan2!( indv, nhtolm, lpl, lpx, ap, qradG, ih, jh, isp, [0.0], ylmk0, qgm )
            qq_nt[ih,jh,isp] = pw.CellVolume * real(qgm[1])
            qq_nt[jh,ih,isp] = pw.CellVolume * real(qgm[1])
        end
    end

    # finally we set the atomic specific qq_at matrices
    for ia in 1:Natoms
        @views qq_at[:,:,ia] = qq_nt[:,:,atm2species[ia]]
    end

    return qradG, qq_nt, qq_at
end




# From aainit of QE-6.6
# lli = lmaxkb + 1
function _calc_clebsch_gordan( lmaxkb::Int64 )

    lli = lmaxkb + 1

    # Those are parameters (HARDCODED)
    lmaxx = 3         # max non local angular momentum (l=0 to lmaxx)
    lqmax = 2*lmaxx + 1  # max number of angular momenta of Q

    # maximum number of combined angular momentum
    nlx = (lmaxx + 1)^2
    # maximum magnetic angular momentum of Q
    mx = 2*lqmax - 1

    llx = (2*lli - 1)^2

    @assert (2*lli-1) <= lqmax
    @assert lli >= 0
    @assert (lli*lli) <= nlx

    r = zeros(Float64, 3, llx)
    for ir in 1:llx
        c = 2.0*rand() - 1.0
        s = sqrt(1.0 - c^2)
        ϕ = 2π * rand()
        r[1,ir] = s * cos(ϕ)
        r[2,ir] = s * sin(ϕ)
        r[3,ir] = c
    end

    # generate the real spherical harmonics for the array: ylm(ir,lm)
    Ylm = zeros(Float64, llx, llx)
     _lmax = round(Int64, sqrt(llx) - 1)
    Ylm_real_qe!(_lmax, r, Ylm) # Ylm_real_qe accept l value starting from 0
    Ylminv = inv(Ylm) # The inverse

    # Clebsch-Gordan coefficients for spherical harmonics
    ap = zeros(Float64, lqmax*lqmax, nlx, nlx)
    # for each pair of combined momenta lm(1),lm(2): 
    lpx = zeros(Int64, nlx, nlx)      # maximum combined angular momentum LM
    lpl = zeros(Int64, nlx, nlx, mx)  # list of combined angular momenta  LM

    # for each li,lj compute ap(l,li,lj) and the indices, lpx and lpl
    #println()
    for li in 1:lli*lli
        #println()
        for lj in 1:lli*lli
            lpx[li,lj] = 0
            for l in 1:llx
                ap[l,li,lj] = _compute_ap(l, li, lj, llx, Ylm, Ylminv)
                if abs(ap[l,li,lj]) > 1e-3
                    lpx[li,lj] = lpx[li,lj] + 1   # increment
                    if lpx[li,lj] > mx
                        println("Error: mx dimension too small: lpx(li,lj)")
                        error("Error in calc_clebsch_gordan")
                    end
                    lpl[li, lj, lpx[li,lj]] = l
                    #@printf("%4d %4d %4d %18.10f\n", l, li, lj, ap[l,li,lj])
                end # if
            end
        end
    end

    return ap, lpx, lpl

end

function _compute_ap(l, li, lj, llx, ylm, mly)
    res = 0.0
    for i in 1:llx
        res += mly[l,i]*ylm[i,li]*ylm[i,lj]
    end
    return res
end


function _calc_qradG(
    pw::PWGrid,
    pspots::Vector{PsPot_UPF}
)
    CellVolume = pw.CellVolume
    Nspecies = length(pspots)
    prefr = 4π/CellVolume
    ecutrho = pw.ecutrho

    ndm = 0
    for isp in 1:Nspecies
        if ndm < pspots[isp].kkbeta
            ndm = pspots[isp].kkbeta
        end
    end

    #qnorm = 0.0 # XXX HARDCODED, no k-points norm of (q + k) ?
    dq = 0.01 # XXX HARDCODED
    cell_factor = 1.0 # hardcoded
    nqxq = floor(Int64, sqrt(2*ecutrho)/dq + 4) # factor of 2 in 2*ecutrho (convert to Ry)

    besr = zeros(Float64, ndm)
    aux = zeros(Float64, ndm)
    # Radial Fourier transform of
    qradG = Array{Array{Float64,3},1}(undef,Nspecies)
    for isp in 1:Nspecies
        Nproj = pspots[isp].Nproj
        Nn2 = round(Int64, Nproj*(Nproj+1)/2)
        # Determine lmaxq
        lmaxkb = -1
        for i in 1:pspots[isp].Nproj
            lmaxkb = max(lmaxkb, pspots[isp].proj_l[i])
        end
        lmaxq = 2*lmaxkb + 1
        qradG[isp] = zeros(Float64, nqxq, Nn2, lmaxq)
    end
    # third dimension can be changed to (lmaxkb[isp] + 1)
    # at least lmaxq=1


    for isp in 1:Nspecies
        #
        psp = pspots[isp] # shorthand
        #
        # skip if not using ultrasoft
        if !pspots[isp].is_ultrasoft
            continue
        end

        for l in 0:psp.nqlc-1
            # note that l is the true (combined) angular momentum
            # and that the arrays have dimensions 0..l (no more 1..l+1)
            for iq in 1:nqxq
                q = (iq - 1) * dq
                # here we compute the spherical bessel function for each q_i
                for ir in 1:psp.kkbeta
                    besr[ir] = sphericalbesselj(l, q*psp.r[ir])
                end
                #
                for nb in 1:psp.Nproj
                    # the Q are symmetric w.r.t indices
                    for mb in nb:psp.Nproj
                        ijv = round(Int64, mb*(mb - 1)/2) + nb
                        lnb = psp.proj_l[nb]
                        lmb = psp.proj_l[mb]
                        # Selection rule
                        cond1 = l >= abs( lnb - lmb )
                        cond2 = l <=      lnb + lmb
                        cond3 = mod(l + lnb + lmb, 2) == 0
                        # XXX: Check qfuncl for this conditions?
                        if cond1 & cond2 & cond3
                            for ir in 1:psp.kkbeta
                                aux[ir] = besr[ir] * psp.qfuncl[ir,ijv,l+1]
                            end
                            # and then we integrate with all the Q functions
                            qradG[isp][iq,ijv,l+1] = PWDFT.integ_simpson( psp.kkbeta, aux, psp.rab )
                        end
                    end
                end # igl
            end # l
        end
        #qradG[isp][:,:,:] = 4*qradG[isp][:,:,:]*prefr
        # Factor of 4 to fix the unit of Deeq and op_S

        @views qradG[isp][:,:,:] = qradG[isp][:,:,:]*prefr
    end
    return qradG
end


# TODO: 
function calc_betaNL_psi(
    ik::Int64,
    pspotNL::PsPotNL_UPF,
    psi::AbstractArray{ComplexF64}
)

    return pspotNL.betaNL[ik]' * psi
end


import Base: print
function print( io::IO, pspotNL::PsPotNL_UPF )
    
    println("------------")
    println("PsPotNL_UPF:")
    println("------------")
    
    println("lmaxx   = ", pspotNL.lmaxx)
    println("lqmax   = ", pspotNL.lqmax)
    println("lmaxkb  = ", pspotNL.lmaxkb)
    println("NbetaNL = ", pspotNL.NbetaNL)

    return
end