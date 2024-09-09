function calc_stress_Ps_nloc!(
    atoms, pw, pspots, pspotNL, electrons, potentials, psiks,
    stress_Ps_nloc
)

    G = pw.gvec.G
    idx_gw2g = pw.gvecw.idx_gw2g
    Ngw = pw.gvecw.Ngw
    Ngwx = maximum(Ngw)
    kfac = ones(Float64, Ngwx) # assume qcutz=0
    Gk = zeros(Float64, Ngwx, 3)
    Gk_length_inv = zeros(Float64, Ngwx)
    #
    k = pw.gvecw.kpoints.k
    Nkpt = pw.gvecw.kpoints.Nkpt
    wk = pw.gvecw.kpoints.wk
    #
    Nstates = electrons.Nstates
    Nspin = electrons.Nspin
    Focc = electrons.Focc

    fill!(stress_Ps_nloc, 0.0)
    for ispin in 1:Nspin, ik in 1:Nkpt
        #
        ikspin = ik + (ispin - 1)*Nkpt
        #
        for igw in 1:Ngw[ik]
            ig = idx_gw2g[ik][igw]
            for i in 1:3
                Gk[igw,i] = G[i,ig] + k[i,ik]
            end
            # NOTE: reversed the usual dimension of Gk[1:3,1:Ngwx] to Gk[1:Ngwx,1:3]
            q = sqrt( Gk[igw,1]^2 + Gk[igw,2]^2 + Gk[igw,3]^2 )
            if q > 1e-8
                Gk_length_inv[igw] = 1.0 / q
            else
                Gk_length_inv[igw] = 0.0
            end
        end

        _stress_Ps_nloc_k!(
            ispin, ik,
            atoms, pw, pspots, pspotNL, electrons,
            Gk, Gk_length_inv,
            psiks, stress_Ps_nloc)

    end


    # Add USPP contribution
    _add_stress_uspp!(atoms, pw, pspots, pspotNL, potentials, stress_Ps_nloc)



    for l in 1:3, m in 1:(l-1)
        stress_Ps_nloc[m,l] = stress_Ps_nloc[l,m]
    end
    # Scale by -1/pw.CellVolume
    stress_Ps_nloc[:,:] .*= (-1/pw.CellVolume)


    return
end


function _stress_Ps_nloc_k!(
    ispin::Int64, ik::Int64,
    atoms, pw, pspots, pspotNL, electrons, 
    Gk, Gk_length_inv,
    psiks, stress_Ps_nloc
)

    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species

    Nstates = electrons.Nstates
    ebands = electrons.ebands

    nh = pspotNL.nh
    nhm = pspotNL.nhm
    Deeq = pspotNL.Deeq
    indv_ijkb0 = pspotNL.indv_ijkb0
    NbetaNL = pspotNL.NbetaNL

    Nkpt = pw.gvecw.kpoints.Nkpt
    wk = pw.gvecw.kpoints.wk
    Ngwk = pw.gvecw.Ngw[ik]

    ikspin = ik + (ispin - 1)*Nkpt
    psi = psiks[ikspin]
    Vnl_KB = pspotNL.betaNL[ik]
    betaNL_psi = Vnl_KB' * psi  # XXX precalculate this

    fac = 2*wk[ik] # XXX: this is different from QE, it does not depend on ist
    # XXX Probably Focc also should be accounted for (also some normalization by Nstates)
    println("fac = ", fac)

    Deff = zeros(Float64, nhm, nhm, Natoms)

    evps = 0.0
    for ist in 1:Nstates
        if abs(fac) < 1e-9
            continue
        end
        _calc_Deff!( ispin, atoms, pspotNL, ebands[ist,ikspin], Deff )
        @printf("ist = %5d eband = %18.10f, fac = %18.10f\n", ist, ebands[ist,ikspin], fac)
        #println("sum Deff = ", sum(Deff))
        #
        ijkb0 = 0
        for isp in 1:Nspecies, ia in 1:Natoms
            #
            if atm2species[ia] != isp
                continue
            end
            #
            psp = pspots[isp]
            #
            for ih in 1:nh[isp]
                ikb = ijkb0 + ih
                evps += fac * Deff[ih,ih,ia] * abs(betaNL_psi[ikb,ist])^2
                #@printf("%5d %5d %18.10f\n", ia, ih, Deff[ih,ih,ia])
                #@printf("Current value of evps = %18.10f\n", evps)
                #
                if psp.is_ultrasoft || (psp.Nproj > 1)
                    # only in the US case there is a contribution for jh != ih
                    # we use here the symmetry in the interchange of ih and jh
                    for jh in (ih + 1):nh[isp]
                        jkb = ijkb0 + jh
                        bibj = real( conj(betaNL_psi[ikb,ist]) * betaNL_psi[jkb,ist] )
                        evps += Deff[ih,jh,ia] * fac * bibj * 2
                        #@printf("Nondiagonal: current value of evps = %18.10f, Deff = %18.10f\n", evps, Deff[ih,jh,ia])
                    end
                end
            end
            ijkb0 += nh[isp]
        end
    end # Nstates
    #    
    println("evps = ", evps)
    #
    for l in 1:3
        stress_Ps_nloc[l,l] -= evps
    end

    #
    # non diagonal contribution - derivative of the bessel function (radial part?)
    #
    dvkb = zeros(ComplexF64, Ngwk, NbetaNL)
  
    _gen_us_dj!(ik, atoms, pw, pspots, pspotNL, dvkb)
    @printf("ik, sum dvkb Bessel = %5d [%18.10f %18.10f]\n", ik, real(sum(dvkb)), imag(sum(dvkb)))

    work2 = zeros(ComplexF64, Ngwk)
    work1 = zeros(ComplexF64, Ngwk)
    #
    for ist in 1:Nstates
        #
        fill!(work2, 0.0 + im*0.0)
        #
        _calc_Deff!( ispin, atoms, pspotNL, ebands[ist,ikspin], Deff )
        #
        ijkb0 = 0
        for isp in 1:Nspecies, ia in 1:Natoms
            #
            if atm2species[ia] != isp
                continue
            end
            #
            psp = pspots[isp]
            uspp_or_multiproj = psp.is_ultrasoft || (psp.Nproj > 1)
            #
            for ih in 1:nh[isp]
                ikb = ijkb0 + ih
                if !uspp_or_multiproj
                    ps = betaNL_psi[ikb,ist] * Deeq[ih,ih,ia,ispin]
                else
                    # in the USPP case there is a contribution also for jh != ih
                    ps = 0.0 + im*0.0
                    for jh in 1:nh[isp]
                        jkb = ijkb0 + jh
                        ps += betaNL_psi[jkb,ist] * Deff[ih,jh,ia]
                    end
                end
                #
                for igw in 1:Ngwk
                    work2[igw] = ps * dvkb[igw,ikb] + work2[igw]
                end
            end
            ijkb0 = ijkb0 + nh[isp]
        end
        #
        for ipol in 1:3, jpol in 1:ipol
            for igw in 1:Ngwk
                work1[igw] = psi[igw,ist] * Gk[igw,ipol] * Gk[igw,jpol] * Gk_length_inv[igw]
            end
            dd = real(dot(work1, work2))
            stress_Ps_nloc[ipol,jpol] -= 2 * (2 * wk[ik]) * dd # factor of two for spin degeneracy
            # wk should be ist dependent
        end
    end # ist

    println("\nstress_Ps_nloc after deriv Bessel contrib (in Ry/bohr^3) = ")
    for i in 1:3
        @printf("%18.10f %18.10f %18.10f\n", 2*stress_Ps_nloc[i,1], 2*stress_Ps_nloc[i,2], 2*stress_Ps_nloc[i,3])
    end


    #
    # non diagonal contribution - derivative of the spherical harmonics
    # (no contribution from l=0)

    # Immediate return
    if pspotNL.lmaxkb == 0
        return
    end

    xyz = zeros(Float64, 3, 3)
    xyz[1,:] = [1.0, 0.0, 0.0]
    xyz[2,:] = [0.0, 1.0, 0.0]
    xyz[3,:] = [0.0, 0.0, 1.0]
    # also can use Matrix(I(3))

    for ipol in 1:3

        _gen_us_dy!(ik, xyz[:,ipol], atoms, pw, pspots, pspotNL, dvkb)
        
        for ist in  1:Nstates
            #
            fill!(work2, 0.0 + im*0.0)
            _calc_Deff!( ispin, atoms, pspotNL, ebands[ist,ikspin], Deff )

            ijkb0 = 0
            for isp in 1:Nspecies, ia in 1:Natoms
                if atm2species[ia] != isp
                    continue
                end
                psp = pspots[isp]
                uspp_or_multiproj = psp.is_ultrasoft || (psp.Nproj > 1)
                for ih in 1:nh[isp]
                    ikb = ijkb0 + ih
                    if !uspp_or_multiproj
                        ps = betaNL_psi[ikb,ist] * Deeq[ih,ih,ia,ispin]
                    else
                        # in the US case there is a contribution also for jh<>ih
                        ps = 0.0 + im*0.0
                        for jh in 1:nh[isp]
                            jkb = ijkb0 + jh
                            ps += betaNL_psi[jkb,ist] * Deff[ih,jh,ia]
                        end
                    end
                    #
                    for igw in 1:Ngwk
                        work2[igw] += ps * dvkb[igw,ikb]
                    end
                end
                ijkb0 += nh[isp]
            end
            for jpol in 1:ipol
                for igw in 1:Ngwk
                    work1[igw] = psi[igw,ist] * Gk[igw,jpol]
                end
                dd = real(dot(work1, work2))
                stress_Ps_nloc[ipol,jpol] -= 2 * (2 * wk[ik]) * dd
            end
        end # ist
    end # ipol
    println("\nstress_Ps_nloc after deriv YLm contrib (in Ry/bohr^3) = ")
    for i in 1:3
        @printf("%18.10f %18.10f %18.10f\n", 2*stress_Ps_nloc[i,1], 2*stress_Ps_nloc[i,2], 2*stress_Ps_nloc[i,3])
    end
    #
    return

end


function _add_stress_uspp!(atoms, pw::PWGrid, pspots, pspotNL, potentials, stress_Ps_nloc)

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
    
    # Fourier transform of the total effective potential
    Vg = zeros(ComplexF64, Ng, Nspin)
    aux = zeros(ComplexF64, Npoints)
    planfw = plan_fft!( zeros(ComplexF64,pw.Ns) ) # using default plan
    for ispin in 1:Nspin
        aux[:] .= potentials.Total[:,ispin]
        #R_to_G!(pw, aux) # This cannot be used because FFT plan is not read from serialization data
        ff = reshape(aux, pw.Ns)
        planfw*ff
        @views aux[:] /= Npoints # rescale
        #
        for ig in 1:Ng
            ip = idx_g2r[ig]
            Vg[ig,ispin] = aux[ip]
        end
    end

    aux1 = zeros(ComplexF64, Ng, 3)
    aux2 = zeros(ComplexF64, Ng, Nspin)

    lmaxq = 2*lmaxkb + 1 # using 1-indexing
    _lmax = lmaxq - 1 # or 2*lmaxkb
    ylmk0 = zeros(Float64, Ng, lmaxq*lmaxq)
    Ylm_real_qe!(_lmax, G, ylmk0) # Ylm_real_qe accept l value starting from 0

    dylmk0 = zeros(Float64, Ng, lmaxq*lmaxq)

    fac = zeros(Float64, 3, Nspin)
    s_us = zeros(Float64, 3, 3)
    
    # here we compute the integral Q*V for each atom,
    #       I = sum_G i G_a exp(-iR.G) Q_nm v^*
    # (no contribution from G=0)
    #
    for ipol in 1:3
        dYlm_real_qe!(_lmax, G, dylmk0, ipol)
        for isp in 1:Nspecies
            #
            psp = pspots[isp]
            need_augmentation = psp.is_ultrasoft || psp.is_paw
            if !need_augmentation
                continue # skip for this species
            end
            #
            nij = Int64(nh[isp]*(nh[isp] + 1)/2)
            dQgm = zeros(ComplexF64, Ng, nij)
            tbecsum = zeros(Float64, nij, Nspin)
            ijh = 0
            for ih in 1:nh[isp], jh in ih:nh[isp]
                ijh = ijh + 1
                @views dqvan2!( pspotNL, ipol, ih, jh, isp, G, G2, ylmk0, dylmk0, dQgm[:,ijh] )
            end
            #
            for ia in 1:Natoms
                if atm2species[ia] != isp
                    continue
                end
                @views tbecsum[:,:] .= becsum[1:nij,ia,1:Nspin]
                #CALL dgemm( 'N', 'N', 2*ngm_l, nspin, nij, 1.0_dp, &
                # qgm, 2*ngm_l, tbecsum, nij, 0.0_dp, aux2, 2*ngm_l )
                aux2[:,:] = dQgm * tbecsum
                
                for ispin in 1:Nspin, ig in 1:Ng
                    aux2[ig,ispin] *= conj(Vg[ig,ispin])
                end          
                #
                for ig in 1:Ng
                    GX = G[1,ig]*atpos[1,ia] + G[2,ig]*atpos[2,ia] + G[3,ig]*atpos[3,ia]
                    Sf = cos(GX) + im*sin(GX)
                    aux1[ig,1] = Sf * G[1,ig]
                    aux1[ig,2] = Sf * G[2,ig]
                    aux1[ig,3] = Sf * G[3,ig]
                end
                #CALL DGEMM('T','N', 3, nspin, 2*ngm_l, 1.0_dp, aux1, 2*ngm_l, &
                #       aux2, 2*ngm_l, 0.0_dp, fac, 3 )
                fac[:,:] = real(aux1' * aux2)
                for ispin in 1:Nspin, jpol in 1:3
                    s_us[ipol,jpol] -= pw.CellVolume * fac[jpol,ispin]
                end
            end
        end
    end
    # Add
    stress_Ps_nloc[:,:] .+= s_us[:,:]
    # factor of 2 if using gamma only
    #
    return
end
