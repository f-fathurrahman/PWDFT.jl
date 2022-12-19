function calc_rhoe_uspp( Ham::Hamiltonian, psiks::BlochWavefunc )
    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Rhoe = zeros(Float64, Npoints, Nspin)
    calc_rhoe_uspp!(Ham, psiks, Rhoe)
    return Rhoe
end

function calc_rhoe_uspp!(
    Ham::Hamiltonian,
    psiks::BlochWavefunc,
    Rhoe::Array{Float64,2}
)
    pw = Ham.pw
    Focc = Ham.electrons.Focc
    Nspin = Ham.electrons.Nspin
    Nelectrons_true = Ham.electrons.Nelectrons

    CellVolume  = pw.CellVolume
    Ns = pw.Ns
    Nkpt = pw.gvecw.kpoints.Nkpt
    Ngw = pw.gvecw.Ngw
    wk = pw.gvecw.kpoints.wk
    idx_gw2r = pw.gvecw.idx_gw2r
    Npoints = prod(Ns)
    if pw.using_dual_grid
        NpointsSmooth = prod(pw.Nss)
    else
        NpointsSmooth = Npoints
    end
    Nstates = size(psiks[1])[2]

    RhoeSmooth = zeros(Float64, NpointsSmooth, Nspin)
    ctmp = zeros(ComplexF64, NpointsSmooth)

    # dont forget to zero out the Rhoe first
    fill!(Rhoe, 0.0)

    NptsPerSqrtVol = NpointsSmooth/sqrt(CellVolume)

    dVol = CellVolume/Npoints # dense

    nhm = Ham.pspotNL.nhm
    Natoms = Ham.atoms.Natoms
    Nbecsum = Int64( nhm * (nhm + 1)/2 )

    becsum = zeros(Float64, Nbecsum, Natoms, Nspin)

    #
    for ispin in 1:Nspin, ik in 1:Nkpt
        #
        ikspin = ik + (ispin - 1)*Nkpt
        psi = psiks[ikspin]
        #
        for ist in 1:Nstates
            #
            fill!(ctmp, 0.0 + im*0.0)
            #
            for igw in 1:Ngw[ik]
                ip = idx_gw2r[ik][igw]
                ctmp[ip] = psi[igw,ist]
            end
            # to real space (use smooth grid if using dual grid)
            G_to_R!(pw, ctmp, smooth=pw.using_dual_grid)
            # Renormalize
            for ip in 1:NpointsSmooth
                ctmp[ip] *= NptsPerSqrtVol
            end
            #
            w = wk[ik]*Focc[ist,ikspin]
            #
            for ip in 1:NpointsSmooth
                # accumulate
                RhoeSmooth[ip,ispin] += w*real( conj(ctmp[ip])*ctmp[ip] )
            end
        end
    end # ik, ispin

    #println("calc_rhoe_uspp: integ RhoeSmooth = ", sum(RhoeSmooth)*dVol)

    # Interpolate
    for ispin in 1:Nspin
        @views smooth_to_dense!(pw, RhoeSmooth, Rhoe)
    end

    #println("calc_rhoe_uspp: integ Rhoe = ", sum(Rhoe)*dVol)

    #
    # Add ultrasoft contrib
    #
    if pw.using_dual_grid
        fill!(becsum, 0.0) # zero out becsum
        for ispin in 1:Nspin, ik in 1:Nkpt
            _add_becsum!(ik, ispin, Ham, psiks, becsum)
        end
        _add_usdens!(Ham, becsum, Rhoe) # using real space
    end

    #println("calc_rhoe_uspp: integ Rhoe after adding becsum = ", sum(Rhoe)*dVol)

    # renormalize
    #if renormalize
    #    integ_rho = sum(Rhoe)*CellVolume/Npoints
    #    for i in 1:length(Rhoe)
    #        Rhoe[i] *= Nelectrons_true/integ_rho
    #    end
    #end

    # FIXME: Need to calculate RhoeG here?

    # Symmetrize Rhoe if needed
    if Ham.sym_info.Nsyms > 1
        symmetrize_rhoe!( Ham.pw, Ham.sym_info, Ham.rhoe_symmetrizer, Rhoe )
    end

    return
end


function _add_usdens!( Ham, becsum, Rhoe )

    G = Ham.pw.gvec.G
    G2 = Ham.pw.gvec.G2
    Ng = Ham.pw.gvec.Ng
    Nspin = Ham.electrons.Nspin
    lmaxkb = Ham.pspotNL.lmaxkb
    idx_g2r = Ham.pw.gvec.idx_g2r
    nh = Ham.pspotNL.nh

    Nspecies = Ham.atoms.Nspecies
    Natoms = Ham.atoms.Natoms
    atm2species = Ham.atoms.atm2species
    atpos = Ham.atoms.positions

    lmaxq = 2*lmaxkb + 1 # using 1-indexing
    ylmk0 = zeros(Float64, Ng, lmaxq*lmaxq)
    _lmax = lmaxq - 1 # or 2*lmaxkb
    # Ylm_real_qe accept l value starting from 0
    Ylm_real_qe!(_lmax, G, ylmk0)

    aux = zeros(ComplexF64, Ng, Nspin)
    Qgm = zeros(ComplexF64, Ng)

    for isp in 1:Nspecies

        if !Ham.pspots[isp].is_ultrasoft
            continue
        end

        # nij = max number of (ih,jh) pairs per atom type nt
        nij = Int64( nh[isp]*(nh[isp] + 1)/2 )

        # count max number of atoms of type isp
        nab = sum(atm2species .== isp)


        tbecsum = zeros(Float64, nij, nab, Nspin)
        Skk = zeros(ComplexF64, Ng, nab)
        #
        nb = 0
        for ia in 1:Natoms
            if atm2species[ia] != isp
                continue
            end
            nb = nb + 1
            tbecsum[:,nb,:] = becsum[1:nij,ia,1:Nspin]
            for ig in 1:Ng
                GX = atpos[1,ia]*G[1,ig] + atpos[2,ia]*G[2,ig] + atpos[3,ia]*G[3,ig]
                Skk[ig,nb] = cos(GX) - im*sin(GX)
            end
        end


        for ispin in 1:Nspin
            # sum over atoms
            aux2 = Skk * tbecsum[:,:,ispin]'
            # sum over lm indices of Q_{lm}
            ijh = 0
            for ih in 1:nh[isp], jh in ih:nh[isp]
                ijh = ijh + 1
                # qgm is complex here
                qvan2!( Ham.pspotNL, ih, jh, isp, G2, ylmk0, Qgm )
                for ig in 1:Ng
                    aux[ig,ispin] += aux2[ig,ijh]*Qgm[ig]
                end
            end
        end
    end

    Npoints = prod(Ham.pw.Ns)
    CellVolume = Ham.pw.CellVolume
    dVol = CellVolume/Npoints

    ctmp = zeros(ComplexF64, Npoints)
    for ig in 1:Ng
        ip = idx_g2r[ig]
        ctmp[ip] = aux[ig,1]
    end
    # Use the dense grid
    G_to_R!(Ham.pw, ctmp)
    ctmp[:] *= Npoints # rescale
    Rhoe[:,1] += real(ctmp)

    return
end


# Based on compute_becsum and sum_bec in QE
#
# Compute "becsum" = \sum_i w_i <psi_i|beta_l><beta_m|\psi_i> term.
# Output in module uspp and (PAW only) in rho%bec (symmetrized)
# if iflag = 1, weights w_k are re-computed.
#
# or:
#
# \[ \sum_i \langle\psi_i|\beta_l\rangle w_i \langle\beta_m|\psi_i\rangle \]
#
# for point "ik" and, for LSDA, spin "current_spin".  
function _add_becsum!( ik, ispin, Ham, psiks, becsum )

    Natoms = Ham.atoms.Natoms
    Nspecies = Ham.atoms.Nspecies
    atm2species = Ham.atoms.atm2species

    Nspin = Ham.electrons.Nspin
    Focc = Ham.electrons.Focc
    
    nhm = Ham.pspotNL.nhm
    nh = Ham.pspotNL.nh
    indv_ijkb0 = Ham.pspotNL.indv_ijkb0

    Nkpt = Ham.pw.gvecw.kpoints.Nkpt

    ikspin = ik + (ispin-1)*Nkpt
    psi = psiks[ikspin]

    wk = Ham.pw.gvecw.kpoints.wk

    Nstates = size(psi, 2)

    betaNL_psi = Ham.pspotNL.betaNL[ik]' * psi

    for isp in 1:Nspecies

        # Skip if this species does not use ultrasoft
        if !Ham.pspots[isp].is_ultrasoft
            continue
        end
        
        # These are used for GEMM.
        # They can be removed
        auxk1 = zeros(ComplexF64, Nstates, nh[isp])
        auxk2 = zeros(ComplexF64, Nstates, nh[isp])
        #
        # In becp=<vkb_i|psi_j> terms corresponding to atom ia of type isp
        # run from index i=indv_ijkb0[ia]+1 to i=indv_ijkb0[ia] + nh[isp]
        for ia in 1:Natoms
            
            if atm2species[ia] != isp
                continue
            end

            # sum over bands: \sum_i <psi_i|beta_l><beta_m|psi_i> w_i
            # copy into aux1, aux2 the needed data to perform a GEMM
            for ih in 1:nh[isp]
                ikb = indv_ijkb0[ia] + ih
                for ist in 1:Nstates
                    auxk1[ist,ih] = betaNL_psi[ikb,ist]
                    auxk2[ist,ih] = wk[ik]*Focc[ist,ikspin]*betaNL_psi[ikb,ist]
                end
            end
            #
            # only the real part is computed
            #
            #CALL DGEMM ( 'C', 'N', nh(np), nh(np), 2*this_bgrp_nbnd, &
            #     1.0_dp, auxk1, 2*this_bgrp_nbnd, auxk2, 2*this_bgrp_nbnd, &
            #     0.0_dp, aux_gk, nh(np) )
            aux_gk = real(auxk1' * auxk2)

            # TODO: calculated aux_gk directly without auxk1 and auxk2
            
            # copy output from GEMM into desired format
            #
            ijh = 0
            for ih in 1:nh[isp], jh in ih:nh[isp]
                ijh = ijh + 1
                # nondiagonal terms summed and collapsed into a
                # single index (matrix is symmetric wrt (ih,jh))
                if jh == ih
                    becsum[ijh,ia,ispin] += aux_gk[ih,jh]
                else
                    becsum[ijh,ia,ispin] += aux_gk[ih,jh]*2.0
                end
            end
        end
    end
    return 
end
