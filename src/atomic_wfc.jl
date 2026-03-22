# This probably should go to PsPot_UPF
function init_tab_at(psp::PsPot_UPF, pw::PWGrid)

    ecutwfc = pw.ecutwfc
    CellVolume = pw.CellVolume
    cell_factor = 1.0 # XXX HARDCODED
    dq = 0.01 # XXX HARDCODED
    ndm = psp.Nr_rcut
    aux = zeros(Float64, ndm)
    vchi = zeros(Float64, ndm)
    
    # chiq = radial fourier transform of atomic orbitals chi
    pref = 4π / sqrt(CellVolume)
    # needed to normalize atomic wfcs (not a bad idea in general and 
    # necessary to compute correctly lda+U projections)
    
    Nq = floor( Int64, (sqrt(2*ecutwfc)/dq + 4)*cell_factor )
    Nwfc = psp.Nchi
    tab_at = zeros(Float64, Nq, Nwfc)
    for iwf in 1:Nwfc
        if psp.occ_chi[iwf] >= 0.0
            l = psp.lchi[iwf]
            for iq in 1:Nq
                q = dq * (iq - 1)
                PWDFT.qe_sph_bes!(l, q, psp.r[1:ndm], aux)
                for ir in 1:ndm
                    vchi[ir] = psp.chi[ir,iwf] * aux[ir] * psp.r[ir]
                end
                vqint = PWDFT.integ_simpson( ndm, vchi, psp.rab )
                tab_at[iq,iwf] = vqint * pref
            end
        end # if
    end # iwf
    return tab_at
end


function atomic_wfc!( ik, atoms, pspots, pw, wfcatom )

    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species
    Natoms = atoms.Natoms
    atpos = atoms.positions

    Ngwk = pw.gvecw.Ngw[ik]
    idx_gw2g = pw.gvecw.idx_gw2g
    G = pw.gvec.G
    k = pw.gvecw.kpoints.k

    dq = 0.01 # XXX HARDCODED

    # calculate max angular momentum required in wavefunctions
    lchi_max = 0
    Nchi_max = 0
    for isp in 1:Nspecies
        psp = pspots[isp]
        Nchi_max = max(Nchi_max, psp.Nchi)
        for iwf in 1:psp.Nchi
            lchi_max = max( lchi_max, psp.lchi[iwf] )
        end
    end
    
    Gk = zeros(Float64, 3, Ngwk)
    Gk2 = zeros(Float64, Ngwk)
    for igk in 1:Ngwk
        ig = idx_gw2g[ik][igk] # index of Gk in G
        Gk[1,igk] = G[1,ig] + k[1,ik]
        Gk[2,igk] = G[2,ig] + k[2,ik]
        Gk[3,igk] = G[3,ig] + k[3,ik]
        Gk2[igk] = Gk[1,igk]^2 +  Gk[2,igk]^2 + Gk[3,igk]^2
    end

    # Ylm_real_qe accept l value starting from 0 (the actual 'physics' angular momentum number)
    ylm = zeros(Float64, Ngwk, (lchi_max+1)^2)
    Ylm_real_qe!(lchi_max, Gk, ylm)
    
    # chiq = radial fourier transform of atomic orbitals chi
    chiq = zeros(Float64, Ngwk, Nchi_max, Nspecies)
    for isp in 1:Nspecies
        psp = pspots[isp]
        tab_at = init_tab_at(psp, pw)
        for nb in 1:psp.Nchi
            if psp.occ_chi[nb] >= 0.0 # why equal to 0 is also included?
                for igw in 1:Ngwk
                    Gm = sqrt(Gk2[igw])
                    px = Gm/dq - floor(Int64, Gm/dq)
                    ux = 1.0 - px
                    vx = 2.0 - px
                    wx = 3.0 - px
                    i0 = floor(Int64, Gm/dq) + 1
                    i1 = i0 + 1
                    i2 = i0 + 2
                    i3 = i0 + 3
                    chiq[igw,nb,isp] = tab_at[i0,nb] * ux * vx * wx / 6.0 +
                                       tab_at[i1,nb] * px * vx * wx / 2.0 -
                                       tab_at[i2,nb] * px * ux * wx / 2.0 +
                                       tab_at[i3,nb] * px * ux * vx / 6.0
                end # for
            end # if
        end # for
    end

    Sk = zeros(ComplexF64, Ngwk)
    Natomwfc = calc_Natomwfc(atoms, pspots)

    Nstates = size(wfcatom, 2)
    fill!(wfcatom, 0.0 + im*0.0)
    #
    n_starting_wfc = 0
    for ia in 1:Natoms
        for igk in 1:Ngwk
            # XXX use dot product here?
            GkX = atpos[1,ia]*Gk[1,igk] + atpos[2,ia]*Gk[2,igk] + atpos[3,ia]*Gk[3,igk]
            Sk[igk] = cos(GkX) - im*sin(GkX)
        end
        isp = atm2species[ia]
        psp = pspots[isp]
        for nb in 1:psp.Nchi
            if psp.occ_chi[nb] >= 0.0
                l = psp.lchi[nb]
                lphase = im^l
                #  the factor i^l MUST BE PRESENT in order to produce
                #  wavefunctions for k=0 that are real in real space
                #
                # FIXME: Noncollinear case is not yet treated
                for m in 1:(2*l + 1)
                    lm = l^2 + m
                    n_starting_wfc += 1
                    if n_starting_wfc > Natomwfc
                        error("Too many atomic wfcs")
                    end
                    for igw in 1:Ngwk
                        wfcatom[igw,n_starting_wfc] = lphase * Sk[igw] * ylm[igw,lm] * chiq[igw,nb,isp]
                    end
                end # for
            end # if
        end # for nb
    end # for atoms

    #@infiltrate

    #
    if Natomwfc < Nstates
        for ist in (Natomwfc+1):Nstates
            for igw in 1:Ngwk
                z = rand() # random modulus
                θ = 2*pi*rand() # random angle
                num = z*(cos(θ) + im*sin(θ))
                denum = Gk2[igw] + 1.0
                wfcatom[igw,ist] = num/denum
            end
        end
    end

    return
end


function initwfc(Ham)
    psiks = zeros_BlochWavefunc(Ham)
    initwfc!(Ham, psiks)
    return psiks
end

function initwfc!(Ham, psiks; Haux=nothing)

    # prepare atomic_wfc
    Nstates = Ham.electrons.Nstates
    Nspin = Ham.electrons.Nspin_wf
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Natomwfc = calc_Natomwfc(Ham.atoms, Ham.pspots)
    Nwfc_start = max(Natomwfc, Nstates)
    
    Ngw = Ham.pw.gvecw.Ngw
    Ngwx = maximum(Ngw)
    wfcatom = zeros(ComplexF64, Ngwx, Nwfc_start)
    Hsub = zeros(ComplexF64, Nwfc_start, Nwfc_start)
    λ = zeros(Float64, Nwfc_start)
    v = zeros(ComplexF64, Nwfc_start, Nwfc_start)

    #ispin = 1
    #ik = 1
    for ispin in 1:Nspin, ik in 1:Nkpt
        #
        ikspin = ik + (ispin-1)*Nkpt
        # Set this because we will do multiplication with Hamiltonian
        Ham.ik = ik
        Ham.ispin = ispin
        #
        Ngwk = Ngw[ik]
        fill!(wfcatom, 0.0)
        @views atomic_wfc!(ik, Ham.atoms, Ham.pspots, Ham.pw, wfcatom[1:Ngwk,:])
        
        #XXX No need for Hsub?
        #psiks[ikspin][:,:] = wfcatom[1:Ngwk,1:Nstates]

        #@infiltrate

        # Need to orthonormalize?
        @views ortho_sqrt!(Ham, wfcatom[1:Ngwk,:])
        #
        Hsub[:,:] = wfcatom[1:Ngwk,:]' * op_H(Ham, wfcatom[1:Ngwk,:])
        λ[:], v[:,:] = eigen(Hermitian(Hsub))
        
        if !isnothing(Haux)
            for ist in 1:Nstates
                Haux[ikspin][ist,ist] = λ[ist] + randn()
            end
        end
        psiks[ikspin][:,:] = wfcatom[1:Ngwk,1:Nstates]*v[1:Nstates,1:Nstates]
        #
        # Force orthonormalization
        #ortho_sqrt!(Ham, psiks[ikspin][1:Ngwk,:])
    end

    return

end