function calc_forces_Ps_nloc_finite_diff(
    atoms::Atoms,    
    pw::PWGrid,
    pspots::Array{PsPot_GTH},
    electrons::Electrons,
    psiks::BlochWavefunc
)

    Δ = 0.001

    pos_orig = copy(atoms.positions)
    Natoms = atoms.Natoms

    kpoints = pw.gvecw.kpoints

    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume
    dVol = CellVolume/Npoints

    F_Ps_nloc = zeros(3,Natoms)

    for ia = 1:Natoms
    for idir = 1:3
        
        # set to original positions
        atoms.positions[:,:] = pos_orig[:,:]
        
        atoms.positions[idir,ia] = pos_orig[idir,ia] + 0.5*Δ        
        
        pspotNL_plus = PsPotNL( atoms, pw, pspots, kpoints )
        Eplus = calc_E_Ps_nloc( atoms, pw, pspots, electrons, pspotNL_plus, psiks)
        
        atoms.positions[idir,ia] = pos_orig[idir,ia] - 0.5*Δ
        
        pspotNL_minus = PsPotNL( atoms, pw, pspots, kpoints )
        Eminus = calc_E_Ps_nloc( atoms, pw, pspots, electrons, pspotNL_minus, psiks)        

        F_Ps_nloc[idir,ia] = -(Eplus - Eminus)/Δ
    end
    end

    atoms.positions[:,:] = pos_orig[:,:]

    return F_Ps_nloc

end


function calc_E_Ps_nloc(
    atoms::Atoms,
    pw::PWGrid,
    pspots::Array{PsPot_GTH,1},
    electrons::Electrons,
    pspotNL::PsPotNL,
    psiks::BlochWavefunc
)

    Nstates = electrons.Nstates
    Focc = electrons.Focc
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    Pspots = pspots
    prj2beta = pspotNL.prj2beta
    Nkpt = pw.gvecw.kpoints.Nkpt
    wk = pw.gvecw.kpoints.wk
    NbetaNL = pspotNL.NbetaNL
    Nspin = electrons.Nspin

    # calculate E_NL
    E_Ps_nloc = 0.0

    betaNL_psi = zeros(ComplexF64,Nstates,NbetaNL)
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        ikspin = ik + (ispin - 1)*Nkpt
        psi = psiks[ikspin]
        betaNL_psi = calc_betaNL_psi( ik, pspotNL.betaNL, psi )
        for ist = 1:Nstates
            enl1 = 0.0
            for ia = 1:Natoms
                isp = atm2species[ia]
                psp = Pspots[isp]
                for l = 0:psp.lmax
                for m = -l:l
                for iprj = 1:psp.Nproj_l[l+1]
                for jprj = 1:psp.Nproj_l[l+1]
                    ibeta = prj2beta[iprj,ia,l+1,m+psp.lmax+1]
                    jbeta = prj2beta[jprj,ia,l+1,m+psp.lmax+1]
                    hij = psp.h[l+1,iprj,jprj]
                    enl1 = enl1 + hij*real(conj(betaNL_psi[ist,ibeta])*betaNL_psi[ist,jbeta])
                end
                end
                end # m
                end # l
            end
            E_Ps_nloc = E_Ps_nloc + wk[ik]*Focc[ist,ikspin]*enl1
        end
    end
    end

    return E_Ps_nloc

end