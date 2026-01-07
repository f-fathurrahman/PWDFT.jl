function op_V_Ps_nloc( Ham::Hamiltonian, psiks::BlochWavefunc )
    Nstates = size(psiks[1],2) # Nstates should be similar for all Bloch states
    Nspin = Ham.electrons.Nspin_channel
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    out = zeros_BlochWavefunc(Ham)
    
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        out[ikspin] = op_V_Ps_nloc( Ham, psiks[ikspin] )
    end
    return out
end

function op_V_Ps_nloc( Ham::Hamiltonian, psi::AbstractArray{ComplexF64} )
    Vpsi = zeros(ComplexF64, size(psi))
    op_V_Ps_nloc!(Ham, psi, Vpsi)
    return Vpsi
end

function op_V_Ps_nloc!(
    Ham::Hamiltonian,
    psiks::BlochWavefunc,
    Hpsiks::BlochWavefunc
)
    Nstates = size(psiks[1],2)
    Nspin = Ham.electrons.Nspin_channel
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt    
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        i = ik + (ispin - 1)*Nkpt
        op_V_Ps_nloc!( Ham, psiks[i], Hpsiks[i] )
    end
    return
end

function op_V_Ps_nloc!(
    Ham::Hamiltonian{PsPot_GTH},
    psi::AbstractArray{ComplexF64},
    Hpsi::AbstractArray{ComplexF64}
)
    #
    ik = Ham.ik
    # Take `Nstates` to be the size of psi and not from `Ham.electrons.Nstates`.
    Nstates = size(psi,2)
    # first dimension of psi should be Ngw[ik]
    atoms = Ham.atoms
    atm2species = atoms.atm2species
    Natoms = atoms.Natoms
    Pspots = Ham.pspots
    prj2beta = Ham.pspotNL.prj2beta
    betaNL = Ham.pspotNL.betaNL
    Ngw_ik = Ham.pw.gvecw.Ngw[ik]
    #
    betaNL_psi = calc_betaNL_psi( ik, Ham.pspotNL.betaNL, psi )
    #
    for ia = 1:Natoms
        isp = atm2species[ia]
        psp = Pspots[isp]
        for l in 0:psp.lmax
        for m in -l:l
            for iprj in 1:psp.Nproj_l[l+1]
            for jprj in 1:psp.Nproj_l[l+1]
                ibeta = prj2beta[iprj,ia,l+1,m+psp.lmax+1]
                jbeta = prj2beta[jprj,ia,l+1,m+psp.lmax+1]
                hij = psp.h[l+1,iprj,jprj]
                for ist in 1:Nstates, igw in 1:Ngw_ik
                    Hpsi[igw,ist] = Hpsi[igw,ist] + hij*betaNL[ik][igw,ibeta]*betaNL_psi[ist,jbeta]
                end
            end # iprj
            end # jprj
        end # m
        end # l
    end
    return
end


# For UPF pspot
function op_V_Ps_nloc!(
    Ham::Hamiltonian{PsPot_UPF},
    psi::AbstractArray{ComplexF64},
    Hpsi::AbstractArray{ComplexF64}
)

    if Ham.options.noncollinear
        op_V_Ps_nloc_noncollinear!(Ham, psi, Hpsi)
        return
    end
    
    ik = Ham.ik
    ispin = Ham.ispin
    Nstates = size(psi, 2)

    Nspecies = Ham.atoms.Nspecies
    Natoms = Ham.atoms.Natoms
    atm2species = Ham.atoms.atm2species

    NbetaNL = Ham.pspotNL.NbetaNL
    @views Deeq = Ham.pspotNL.Deeq[:,:,:,ispin]
    nh = Ham.pspotNL.nh
    indv_ijkb0 = Ham.pspotNL.indv_ijkb0
    betaNL_k = Ham.pspotNL.betaNL[ik]

    # Use * directly?
    betaNL_psi = calc_betaNL_psi(ik, Ham.pspotNL, psi)

    ps = zeros(ComplexF64, NbetaNL, Nstates)

    for isp in 1:Nspecies
        #
        if nh[isp] == 0
            continue
        end
        #
        for ia in 1:Natoms
            if atm2species[ia] != isp
                continue
            end
            idx1 = (indv_ijkb0[ia]+1):(indv_ijkb0[ia]+nh[isp])
            @views ps[idx1,:] = Deeq[1:nh[isp],1:nh[isp],ia] * betaNL_psi[idx1,:]
        end
    end

    # Accumulate
    @views Hpsi[:,:] += betaNL_k * ps

    return

end


function op_V_Ps_nloc_noncollinear!(
    Ham::Hamiltonian{PsPot_UPF},
    psi::AbstractArray{ComplexF64},
    Hpsi::AbstractArray{ComplexF64}
)
    ik = Ham.ik
    ispin = Ham.ispin
    @assert ispin == 1
    Npol = 2 # it must be this
    Nstates = size(psi, 2)
    #
    Nspecies = Ham.atoms.Nspecies
    Natoms = Ham.atoms.Natoms
    atm2species = Ham.atoms.atm2species
    #
    NbetaNL = Ham.pspotNL.NbetaNL
    Deeq_nc = Ham.pspotNL.Deeq_nc
    nh = Ham.pspotNL.nh
    indv_ijkb0 = Ham.pspotNL.indv_ijkb0
    betaNL_k = Ham.pspotNL.betaNL[ik]
    Ngw_ik = Ham.pw.gvecw.Ngw[ik]
    #
    # Use * directly?
    betaNL_psi = calc_betaNL_psi_noncollinear(ik, Ham.pspotNL, psi)
    ps = zeros(ComplexF64, NbetaNL, Npol, Nstates)
    # XXX: This can be simplified to be just loop over Natoms?
    for isp in 1:Nspecies
        if nh[isp] == 0
            continue
        end
        for ia in 1:Natoms
            if atm2species[ia] != isp
                continue
            end
            for ist in 1:Nstates
                for jh in 1:nh[isp]
                    jkb = indv_ijkb0[ia] + jh
                    for ih in 1:nh[isp]
                        ikb = indv_ijkb0[ia] + ih
                        ps[ikb,1,ist] += Deeq_nc[ih,jh,ia,1]*betaNL_psi[jkb,1,ist] + 
                                          Deeq_nc[ih,jh,ia,2]*betaNL_psi[jkb,2,ist] 
                        ps[ikb,2,ist] += Deeq_nc[ih,jh,ia,3]*betaNL_psi[jkb,1,ist] +
                                          Deeq_nc[ih,jh,ia,4]*betaNL_psi[jkb,2,ist]
                    end
                end
            end
#=
            idx1 = (indv_ijkb0[ia]+1):(indv_ijkb0[ia]+nh[isp])
            # for spin-orbit
            @views ps[idx1,1,:] = Deeq_nc[1:nh[isp],1:nh[isp],ia,1] * betaNL_psi[idx1,1,:]
            @views ps[idx1,2,:] = Deeq_nc[1:nh[isp],1:nh[isp],ia,4] * betaNL_psi[idx1,2,:]
=#               
        end # loop over Natoms
    end
    Hpsir = reshape(Hpsi, Ngw_ik, Npol, Nstates)
    # Accumulate
    @views Hpsir[:,1,:] .+= betaNL_k * ps[:,1,:]
    @views Hpsir[:,2,:] .+= betaNL_k * ps[:,2,:]
    return
end