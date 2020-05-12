import PWDFT: calc_E_kin, calc_E_local, calc_E_Ps_nloc, calc_energies

function calc_E_Ps_nloc( Ham::HamiltonianGamma, psis::BlochWavefuncGamma )

    Nstates = Ham.electrons.Nstates
    Focc = Ham.electrons.Focc
    Nspin = Ham.electrons.Nspin

    Natoms = Ham.atoms.Natoms
    atm2species = Ham.atoms.atm2species
    pspots = Ham.pspots
    prj2beta = Ham.pspotNL.prj2beta
    NbetaNL = Ham.pspotNL.NbetaNL

    # calculate E_NL
    E_Ps_nloc = 0.0

    betaNL_psi = zeros(ComplexF64,Nstates,NbetaNL)
    
    for ispin in 1:Nspin
        
        psi = psis.data[ispin]    
        betaNL_psi = calc_betaNL_psi( Ham.pspotNL, psi )
        
        for ist in 1:Nstates
            enl1 = 0.0
            for ia = 1:Natoms
                isp = atm2species[ia]
                psp = pspots[isp]
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
            E_Ps_nloc = E_Ps_nloc + Focc[ist,ispin]*enl1
        end
    end

    return E_Ps_nloc

end


function calc_E_local( Ham::HamiltonianGamma )

    Npoints = prod(Ham.pw.Ns)
    dVol = Ham.pw.CellVolume/Npoints
    Nspin = Ham.electrons.Nspin
    potentials = Ham.potentials

    Rhoe_tot = zeros(Npoints)
    for ispin = 1:Nspin
        Rhoe_tot[:] = Rhoe_tot[:] + Ham.rhoe[:,ispin]
    end

    E_Hartree = 0.5*dot( potentials.Hartree, Rhoe_tot ) * dVol
    E_Ps_loc = dot( potentials.Ps_loc, Rhoe_tot ) * dVol

    if Ham.xcfunc == "PBE"
        epsxc = calc_epsxc_PBE( Ham.xc_calc, Ham.pw, Ham.rhoe )
    else
        epsxc = calc_epsxc_VWN( Ham.xc_calc, Ham.rhoe )
    end
    E_xc = dot( epsxc, Rhoe_tot ) * dVol

    return E_Ps_loc, E_Hartree, E_xc
end


function calc_E_kin( Ham::HamiltonianGamma, psis::BlochWavefuncGamma )

    Focc = Ham.electrons.Focc
    Nstates = Ham.electrons.Nstates
    Nspin = Ham.electrons.Nspin

    Ngw = Ham.pw.gvecw.Ngw
    idx_gw2g = Ham.pw.gvecw.idx_gw2g
    G2 = Ham.pw.gvec.G2

    E_kin = 0.0
    for ispin = 1:Nspin
        psi = psis.data[ispin]
        for ist = 1:Nstates
            psiKpsi = 0.0
            for igw = 2:Ngw
                ig = idx_gw2g[igw]
                psiKpsi = psiKpsi + abs(psi[ig,ist])^2*G2[ig]
            end
            E_kin = E_kin + Focc[ist,ispin]*psiKpsi
        end
    end
    return E_kin # no factor 0.5, gamma only version

end


function calc_energies( Ham::HamiltonianGamma, psis::BlochWavefuncGamma )
    
    E_kin = calc_E_kin( Ham, psis )

    E_Ps_loc, E_Hartree, E_xc = calc_E_local( Ham )

    if Ham.pspotNL.NbetaNL > 0
        E_Ps_nloc = calc_E_Ps_nloc( Ham, psis )
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

    return energies
end

