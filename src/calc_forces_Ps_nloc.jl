function calc_forces_Ps_nloc( Ham, psiks )
    F_Ps_nloc = zeros(Float64, 3, Ham.atoms.Natoms)
    calc_forces_Ps_nloc!( Ham, psiks, F_Ps_nloc )
    return F_Ps_nloc
end

# This is for PsPot_GTH, no need to pass Ham.potentials
function calc_forces_Ps_nloc!(
    Ham::Hamiltonian{PsPot_GTH},
    psiks::BlochWavefunc,
    F_Ps_nloc::Array{Float64,2}
)
    calc_forces_Ps_nloc!( Ham.atoms, Ham.pw, Ham.pspots, Ham.electrons, Ham.pspotNL, psiks, F_Ps_nloc)
    return
end

# This is for PsPot_UPF, we need to pass Ham.potentials (related to Deff,
# integral between potentials and projectors)
function calc_forces_Ps_nloc!(
    Ham::Hamiltonian{PsPot_UPF},
    psiks::BlochWavefunc,
    F_Ps_nloc::Array{Float64,2}
)
    calc_forces_Ps_nloc!( Ham.atoms, Ham.pw, Ham.pspots, Ham.electrons, Ham.pspotNL,
        Ham.potentials, psiks, F_Ps_nloc)
    return
end


# The actual implementation for PsPot_UPF is here
function calc_forces_Ps_nloc!(
    atoms::Atoms,    
    pw::PWGrid,
    pspots::Vector{PsPot_UPF},
    electrons::Electrons,
    pspotNL::PsPotNL_UPF,
    potentials::Potentials,
    psiks::BlochWavefunc,
    F_Ps_nloc::Matrix{Float64}
)

    betaNL = pspotNL.betaNL
    idx_gw2g = pw.gvecw.idx_gw2g
    G = pw.gvec.G
    NbetaNL = pspotNL.NbetaNL
    Nstates = electrons.Nstates
    Nspin = electrons.Nspin
    Nkpt = pw.gvecw.kpoints.Nkpt

    Ngw = pw.gvecw.Ngw
    Ngwx = maximum(Ngw)
    # To save memory usage
    dbetaNL = zeros(ComplexF64,Ngwx,NbetaNL)

    betaNL_psi = zeros(ComplexF64,NbetaNL,Nstates)
    dbetaNL_psi = zeros(ComplexF64,NbetaNL,Nstates)

    fill!(F_Ps_nloc, 0.0)

    for ispin in 1:Nspin, ik in 1:Nkpt
        #
        ikspin = ik + (ispin - 1)*Nkpt
        psi = psiks[ikspin]
        Ngw_ik = Ngw[ik]
        #
        betaNL_psi[:,:] .= betaNL[ik]' * psi
        # Calculate derivative of betaNL in G-space
        for ipol in 1:3
            for ibeta in 1:NbetaNL, igw in 1:Ngw_ik
                ig = idx_gw2g[ik][igw]
                dbetaNL[igw,ibeta] = -im * betaNL[ik][igw,ibeta] * G[ipol,ig]
            end
            # betaNL psi 
            @views dbetaNL_psi[:,:] .= dbetaNL[1:Ngw_ik,:]' * psi
            #
            # this will call sum over bands
            _force_Ps_nloc_k!(ipol, ik, ispin, atoms, pw, pspots,
                electrons, pspotNL, betaNL_psi, dbetaNL_psi, F_Ps_nloc)
        end
    end

    # This does not depend on k-points
    _add_F_uspp!(atoms, pw, pspots, pspotNL, potentials, F_Ps_nloc)

    return

end


function _force_Ps_nloc_k!(ipol, ik, ispin, 
    atoms::Atoms,
    pw::PWGrid,
    pspots::Vector{PsPot_UPF},
    electrons::Electrons,
    pspotNL::PsPotNL_UPF,
    betaNL_psi, dbetaNL_psi,
    F_Ps_nloc
)

    atm2species = atoms.atm2species
    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies

    ebands = electrons.ebands
    Nstates = electrons.Nstates
    Focc = electrons.Focc

    nh = pspotNL.nh
    nhm = pspotNL.nhm
    Deeq = pspotNL.Deeq
    indv_ijkb0 = pspotNL.indv_ijkb0

    Nkpt = pw.gvecw.kpoints.Nkpt
    wk = pw.gvecw.kpoints.wk

    Deff = zeros(Float64, nhm, nhm, Natoms)

    ikspin = ik + (ispin - 1)*Nkpt
    #
    for ist in 1:Nstates
        #
        _calc_Deff!( ispin, atoms, pspotNL, ebands[ist,ik], Deff )
        #
        fac = wk[ik] * Focc[ist,ikspin]
        #
        for isp in 1:Nspecies, ia in 1:Natoms
            
            if atm2species[ia] != isp
                continue
            end

            psp = pspots[isp]            
            ijkb0 = indv_ijkb0[ia]

            for ih in 1:nh[isp]
                ikb = ijkb0 + ih
                F_Ps_nloc[ipol,ia] -= 2.0*fac*Deff[ih,ih,ia] *
                   real( conj(dbetaNL_psi[ikb,ist]) * betaNL_psi[ikb,ist] )
            end
            
            uspp_or_multiproj = psp.is_ultrasoft || (psp.Nproj > 1)

            if uspp_or_multiproj
                # this case is almost always true for our case
                for ih in 1:nh[isp]
                    ikb = ijkb0 + ih
                    # in US case there is a contribution for jh /= ih. 
                    # We use here the symmetry in the interchange  of ih and jh
                    for jh in (ih+1):nh[isp]
                        jkb = ijkb0 + jh
                        F_Ps_nloc[ipol,ia] -= 2.0*fac*Deff[ih,jh,ia] * 
                            real( conj(dbetaNL_psi[ikb,ist]) * betaNL_psi[jkb,ist] +
                                  dbetaNL_psi[jkb,ist] * conj(betaNL_psi[ikb,ist]) )
                    end
                end
            end # uspp_or_multiproj
        end
    end
    return
end


# This routine computes the contribution to atomic forces due
# to the dependence of the Q function on the atomic position.
# \[ F_{j,\text{at}} = \sum_G \sum_{lm} iG_j\ \text{exp}(-iG*R_\text{at})
#    V^*(G)\ Q_{lm}(G)\ \text{becsum}(lm,\text{at}) \]
# where:
# \[ \text{becsum}(lm,\text{at}) = \sum_i \langle \psi_i|\beta_l\rangle
#    w_i\langle \beta_m|\psi_i\rangle \]
# On output: the contribution is added to \(\text{forcenl}\).
function _add_F_uspp!(atoms, pw::PWGrid, pspots, pspotNL, potentials, F_Ps_nloc)

    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species
    atpos = atoms.positions

    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r
    Npoints = prod(pw.Ns)
    G = pw.gvec.G
    G2 = pw.gvec.G2

    lmaxkb = pspotNL.lmaxkb
    nh = pspotNL.nh
    becsum = pspotNL.becsum

    Nspin = size(potentials.Total, 2)

    ok_uspp_or_paw = any(pspotNL.are_ultrasoft) || any(pspotNL.are_paw)
    if !ok_uspp_or_paw
        return
    end

    F_uspp = zeros(Float64, 3, Natoms)
    fact = 1.0*pw.CellVolume # this is 2*pw.CellVolume if using gamma only

    # Fourier transform of the total effective potential
    Vg = zeros(ComplexF64, Ng, Nspin)
    aux = zeros(ComplexF64, Npoints)
    for ispin in 1:Nspin
        aux[:] .= potentials.Total[:,ispin]
        R_to_G!(pw, aux) #XXX why this will segfault?
        @views aux[:] /= Npoints # rescale
        # Note the factors -i and 2pi/a *units of G) here in V(G)
        for ig in 1:Ng
            ip = idx_g2r[ig]
            Vg[ig,ispin] = -im*aux[ip] # XXX need factor tpiba?
        end
    end
    # Finished calculation -im*V_eff(G)


    lmaxq = 2*lmaxkb + 1 # using 1-indexing
    ylmk0 = zeros(Float64, Ng, lmaxq*lmaxq)
    _lmax = lmaxq - 1 # or 2*lmaxkb
    Ylm_real_qe!(_lmax, G, ylmk0) # Ylm_real_qe accept l value starting from 0
    
    for isp in 1:Nspecies

        psp = pspots[isp]
        need_augmentation = psp.is_ultrasoft || psp.is_paw
        if !need_augmentation
            continue # skip for this species
        end

        # nij = max number of (ih,jh) pairs per atom type nt
        nij = Int64(nh[isp]*(nh[isp] + 1)/2)
        Qgm = zeros(ComplexF64, Ng, nij)
        ijh = 0
        for ih in 1:nh[isp], jh in ih:nh[isp]
            ijh = ijh + 1
            @views qvan2!( pspotNL, ih, jh, isp, G2, ylmk0, Qgm[:,ijh] )
        end
        #
        # nab = number of atoms of type nt
        nab = sum(atm2species .== isp)
        #
        aux1 = zeros(ComplexF64, Ng, nab, 3)
        dDeeq = zeros(Float64, nij, nab, 3, Nspin)
        #
        for ispin in 1:Nspin
            nb = 0
            for ia in 1:Natoms
                if atm2species[ia] != isp
                    continue
                end
                nb = nb + 1
                # aux1 = product of potential, structure factor and iG
                for ig in 1:Ng
                    GX = G[1,ig]*atpos[1,ia] + G[2,ig]*atpos[2,ia] + G[3,ig]*atpos[3,ia]
                    Sf = cos(GX) + im*sin(GX)
                    cfac = Vg[ig,ispin] * Sf
                    aux1[ig,nb,1] = G[1,ig] * cfac
                    aux1[ig,nb,2] = G[2,ig] * cfac
                    aux1[ig,nb,3] = G[3,ig] * cfac
                end
            end
            # dDeeq = dot product of aux1 with the Q functions
            # No need for special treatment of the G=0 term (is zero)
            for ipol in 1:3
                @views dDeeq[:,:,ipol,ispin] .= fact * real(Qgm' * aux1[:,:,ipol]) # XXX
            end
        end
        #
        for ispin in 1:Nspin
            nb = 0
            for ia in 1:Natoms
                if atm2species[ia] != isp
                    continue # skip
                end
                nb = nb + 1
                for ipol in 1:3, ijh in 1:nij
                    F_uspp[ipol,ia] += dDeeq[ijh,nb,ipol,ispin] * becsum[ijh,ia,ispin]
                end
            end
        end
    end
    # Add F_uspp to the output array
    @views F_Ps_nloc[:,:] .+= F_uspp[:,:]
    return
end







# The actual implementation for PsPot_GTH is here
function calc_forces_Ps_nloc!(
    atoms::Atoms,    
    pw::PWGrid,
    pspots::Vector{PsPot_GTH},
    electrons::Electrons,
    pspotNL::PsPotNL,
    psiks::BlochWavefunc,
    F_Ps_nloc::Array{Float64,2}
)

    Nstates = electrons.Nstates
    Focc = electrons.Focc
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    prj2beta = pspotNL.prj2beta
    Nkpt = pw.gvecw.kpoints.Nkpt
    wk = pw.gvecw.kpoints.wk
    NbetaNL = pspotNL.NbetaNL
    Nspin = electrons.Nspin

    if NbetaNL == 0
        return F_Ps_nloc
    end

    betaNL_psi = zeros(ComplexF64,Nstates,NbetaNL)
    dbetaNL_psi = zeros(ComplexF64,3,Nstates,NbetaNL)

    dbetaNL = calc_dbetaNL(atoms, pw, pspots, pspotNL)

    for ispin in 1:Nspin, ik in 1:Nkpt

        ikspin = ik + (ispin - 1)*Nkpt
        psi = psiks[ikspin]
        betaNL_psi = calc_betaNL_psi( ik, pspotNL.betaNL, psi )
        dbetaNL_psi = calc_dbetaNL_psi( ik, dbetaNL, psi )
        
        for ist in 1:Nstates, ia in 1:Natoms
            isp = atm2species[ia]
            psp = pspots[isp]
            for l in 0:psp.lmax, m in -l:l
                for iprj in 1:psp.Nproj_l[l+1], jprj in 1:psp.Nproj_l[l+1]
                    ibeta = prj2beta[iprj,ia,l+1,m+psp.lmax+1]
                    jbeta = prj2beta[jprj,ia,l+1,m+psp.lmax+1]
                    hij = psp.h[l+1,iprj,jprj]
                    fnl1 = hij*conj(betaNL_psi[ist,ibeta])*dbetaNL_psi[:,ist,jbeta]
                    @views F_Ps_nloc[:,ia] += 2*wk[ik]*Focc[ist,ikspin]*real(fnl1)
                end
            end # l,m
        end
    end

    return
end



function calc_dbetaNL(
    atoms::Atoms,
    pw::PWGrid,
    pspots::Vector{PsPot_GTH},
    pspotNL::PsPotNL
)

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    atpos = atoms.positions
    kpoints = pw.gvecw.kpoints
    Nkpt = kpoints.Nkpt
    k = kpoints.k
    Ngw = pw.gvecw.Ngw
    Ngwx = pw.gvecw.Ngwx
    NbetaNL = pspotNL.NbetaNL
    dbetaNL = zeros(ComplexF64, 3, Ngwx, NbetaNL, Nkpt) # make this similar to betaNL?
    G = pw.gvec.G
    g = zeros(Float64,3)
    #
    for ik in 1:Nkpt
        ibeta = 0
        for ia in 1:Natoms
            isp = atm2species[ia]
            psp = pspots[isp]
            for l in 0:psp.lmax, m in -l:l, iprj in 1:psp.Nproj_l[l+1]
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
                    Ylmg = Ylm_real(l,m,g)
                    Pg = eval_proj_G(psp,l,iprj,Gm,pw.CellVolume)
                    @views dbetaNL[:,igk,ibeta,ik] = (-1.0*im)^l * Ylmg * Pg * Sf * im * g[:]
                end
            end
        end
    end  # kpoints
    return dbetaNL
end

function calc_dbetaNL_psi(
    ik::Int64,
    dbetaNL::Array{ComplexF64,4},
    psi::Array{ComplexF64,2}
)
    Ngw_ik = size(psi,1)
    Nstates = size(psi,2)
    NbetaNL = size(dbetaNL,3)  # NOTE: 3rd dimension
    dbetaNL_psi = zeros(ComplexF64, 3, Nstates, NbetaNL)
    #
    # FIXME: This can be made into matrix multiplication instead
    #
    for ist in 1:Nstates, ibeta in 1:NbetaNL
        @views dbetaNL_psi[1,ist,ibeta] = dot( dbetaNL[1,1:Ngw_ik,ibeta,ik], psi[:,ist] )
        @views dbetaNL_psi[2,ist,ibeta] = dot( dbetaNL[2,1:Ngw_ik,ibeta,ik], psi[:,ist] )
        @views dbetaNL_psi[3,ist,ibeta] = dot( dbetaNL[3,1:Ngw_ik,ibeta,ik], psi[:,ist] )
    end
    return dbetaNL_psi
end