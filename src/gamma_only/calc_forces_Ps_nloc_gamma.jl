function calc_forces_Ps_nloc!(
    Ham::HamiltonianGamma,
    psiks::BlochWavefuncGamma,
    F_Ps_nloc::Array{Float64,2}
)
    calc_forces_Ps_nloc!( Ham.atoms, Ham.pw, Ham.pspots, Ham.electrons, Ham.pspotNL, psiks, F_Ps_nloc)
    return
end

function calc_forces_Ps_nloc!(
    atoms::Atoms,    
    pw::PWGridGamma,
    pspots::Array{PsPot_GTH,1},
    electrons::Electrons,
    pspotNL::PsPotNLGamma,
    psis::BlochWavefuncGamma,
    F_Ps_nloc::Array{Float64,2}
)

    Nstates = electrons.Nstates
    Focc = electrons.Focc
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    prj2beta = pspotNL.prj2beta
    NbetaNL = pspotNL.NbetaNL
    Nspin = electrons.Nspin_channel

    if NbetaNL == 0
        return F_Ps_nloc
    end

    betaNL_psi = zeros(ComplexF64,Nstates,NbetaNL)
    dbetaNL_psi = zeros(ComplexF64,3,Nstates,NbetaNL)

    dbetaNL = calc_dbetaNL(atoms, pw, pspots, pspotNL)

    for ispin in 1:Nspin
        psi = psis.data[ispin]
        betaNL_psi = calc_betaNL_psi(pspotNL, psi)
        dbetaNL_psi = calc_dbetaNL_psi_gamma(dbetaNL, psi)
        for ist in 1:Nstates
            for ia in 1:Natoms
                isp = atm2species[ia]
                psp = pspots[isp]
                for l in 0:psp.lmax, m in -l:l
                    for iprj in 1:psp.Nproj_l[l+1], jprj in 1:psp.Nproj_l[l+1]
                        ibeta = prj2beta[iprj,ia,l+1,m+psp.lmax+1]
                        jbeta = prj2beta[jprj,ia,l+1,m+psp.lmax+1]
                        hij = psp.h[l+1,iprj,jprj]
                        fnl1 = hij*conj(betaNL_psi[ist,ibeta])*dbetaNL_psi[:,ist,jbeta]
                        F_Ps_nloc[:,ia] = F_Ps_nloc[:,ia] + 2*Focc[ist,ispin]*real(fnl1)
                    end
                end # l
            end
        end
    end

    return
end



function calc_dbetaNL(
    atoms::Atoms,
    pw::PWGridGamma,
    pspots::Array{PsPot_GTH,1},
    pspotNL::PsPotNLGamma
)

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    atpos = atoms.positions

    Ngw = pw.gvecw.Ngw
    NbetaNL = pspotNL.NbetaNL

    dbetaNL = zeros(ComplexF64,3,Ngw,NbetaNL)

    G = pw.gvec.G
    G2 = pw.gvec.G2
    g = zeros(3)
    Ω = pw.CellVolume
    
    idx_gw2g = pw.gvecw.idx_gw2g
    ibeta = 0
    for ia in 1:Natoms
        isp = atm2species[ia]
        psp = pspots[isp]
        for l in 0:psp.lmax, m in -l:l
        for iprj in 1:psp.Nproj_l[l+1]
            ibeta = ibeta + 1
            for igk in 1:Ngw
                ig = idx_gw2g[igk]
                g[1] = G[1,ig]
                g[2] = G[2,ig]
                g[3] = G[3,ig]
                Gm = sqrt(G2[ig])
                GX = atpos[1,ia]*g[1] + atpos[2,ia]*g[2] + atpos[3,ia]*g[3]
                Sf = cos(GX) - im*sin(GX)
                dbetaNL[:,igk,ibeta] =
                (-1.0*im)^l * Ylm_real(l,m,g)*eval_proj_G(psp,l,iprj,Gm,Ω)*Sf*im*g[:]
            end
        end
        end
    end

    return dbetaNL
end

function calc_dbetaNL_psi_gamma( dbetaNL::Array{ComplexF64,3}, psi::Array{ComplexF64,2} )
    Nstates = size(psi)[2]
    NbetaNL = size(dbetaNL)[3]  # NOTE: 3rd dimension
    dbetaNL_psi = zeros(ComplexF64,3,Nstates,NbetaNL)
    for ist in 1:Nstates
        for ibeta in 1:NbetaNL
            #
            c = dot(dbetaNL[1,:,ibeta], psi[:,ist])
            v = conj(dbetaNL[1,1,ibeta])*psi[1,ist]
            dbetaNL_psi[1,ist,ibeta] = c + conj(c) - v
            #
            c = dot(dbetaNL[2,:,ibeta], psi[:,ist])
            v = conj(dbetaNL[2,1,ibeta])*psi[1,ist]
            dbetaNL_psi[2,ist,ibeta] = c + conj(c) - v
            #
            c = dot(dbetaNL[3,:,ibeta], psi[:,ist])
            v = conj(dbetaNL[3,1,ibeta])*psi[1,ist]
            dbetaNL_psi[3,ist,ibeta] = c + conj(c) - v
        end
    end
    return dbetaNL_psi
end