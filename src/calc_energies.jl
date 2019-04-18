function calc_E_Ps_nloc( Ham::Hamiltonian, psiks::BlochWavefunc )

    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc
    Natoms = Ham.atoms.Natoms
    atm2species = Ham.atoms.atm2species
    Pspots = Ham.pspots
    prj2beta = Ham.pspotNL.prj2beta
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    wk = Ham.pw.gvecw.kpoints.wk
    NbetaNL = Ham.pspotNL.NbetaNL
    Nspin = Ham.electrons.Nspin

    # calculate E_NL
    E_Ps_nloc = 0.0

    betaNL_psi = zeros(ComplexF64,Nstates,NbetaNL)
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        ikspin = ik + (ispin - 1)*Nkpt
        psi = psiks[ikspin]
        betaNL_psi = calc_betaNL_psi( ik, Ham.pspotNL.betaNL, psi )
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


#
# psi is assumed to be already orthonormalized elsewhere
# `potentials` and `Rhoe` are not updated
# Ham is assumed to be already updated at input psi
#
# Ham.energies.NN should be calculated outside this function
function calc_energies( Ham::Hamiltonian, psiks::BlochWavefunc )

    pw = Ham.pw
    potentials = Ham.potentials
    Focc = Ham.electrons.Focc

    CellVolume = pw.CellVolume
    Ns = pw.Ns
    Npoints = prod(Ns)
    dVol = CellVolume/Npoints
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nstates = Ham.electrons.Nstates
    wk = Ham.pw.gvecw.kpoints.wk
    Nspin = Ham.electrons.Nspin

    Ngw = Ham.pw.gvecw.Ngw
    idx_gw2g = Ham.pw.gvecw.idx_gw2g
    G = Ham.pw.gvec.G
    k = Ham.pw.gvecw.kpoints.k
    
    #
    # Kinetic energy
    #
    E_kin = 0.0
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        psi = psiks[ikspin]
        for ist = 1:Nstates
            psiKpsi = 0.0
            for igk = 1:Ngw[ik]
                ig = idx_gw2g[ik][igk]
                Gw2 = (G[1,ig] + k[1,ik])^2 + (G[2,ig] + k[2,ik])^2 + (G[3,ig] + k[3,ik])^2
                psiKpsi = psiKpsi + abs(psi[igk,ist])^2*Gw2
            end
            E_kin = E_kin + wk[ik]*Focc[ist,ikspin]*psiKpsi
        end
    end
    end
    E_kin = 0.5*E_kin

    Rhoe_tot = zeros(Npoints)
    for ispin = 1:Nspin
        Rhoe_tot[:] = Rhoe_tot[:] + Ham.rhoe[:,ispin]
    end

    cRhoeG = conj(R_to_G(pw, Rhoe_tot))/Npoints
    V_HartreeG = R_to_G(pw, potentials.Hartree)
    V_Ps_locG = R_to_G(pw, potentials.Ps_loc)

    E_Hartree = 0.0
    E_Ps_loc = 0.0
    for ig = 2:pw.gvec.Ng
        ip = pw.gvec.idx_g2r[ig]
        E_Hartree = E_Hartree + abs( V_HartreeG[ip]*cRhoeG[ip] )
        E_Ps_loc = E_Ps_loc + real( V_Ps_locG[ip]*cRhoeG[ip] )
    end
    #E_Ps_loc = E_Ps_loc + real( V_Ps_locG[1]*cRhoeG[1] ) # should be the same as pspcore ene)
    E_Hartree = 0.5*E_Hartree*dVol
    E_Ps_loc = E_Ps_loc*dVol

    if Ham.xcfunc == "PBE"
        epsxc = calc_epsxc_PBE( Ham.pw, Ham.rhoe )
    else
        epsxc = calc_epsxc_VWN( Ham.rhoe )
    end
    E_xc = dot( epsxc, Rhoe_tot ) * dVol

    if Ham.pspotNL.NbetaNL > 0
        E_Ps_nloc = calc_E_Ps_nloc( Ham, psiks )
    else
        E_Ps_nloc = 0.0
    end

    energies = Energies()
    energies.Kinetic = E_kin
    energies.Ps_loc  = E_Ps_loc
    energies.Ps_nloc = E_Ps_nloc
    energies.Hartree = E_Hartree
    energies.XC      = E_xc
    energies.NN      = Ham.energies.NN
    energies.PspCore = Ham.energies.PspCore

    return energies
end

