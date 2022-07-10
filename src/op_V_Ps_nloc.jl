function op_V_Ps_nloc( Ham::Hamiltonian, psiks::BlochWavefunc )
    Nstates = size(psiks[1],2) # Nstates should be similar for all Bloch states
    Nspin = Ham.electrons.Nspin
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
    Nspin = Ham.electrons.Nspin
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
    Ham::Hamiltonian{Txc,PsPot_GTH},
    psi::AbstractArray{ComplexF64},
    Hpsi::AbstractArray{ComplexF64}
) where Txc <: AbstractXCCalculator
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
    Ham::Hamiltonian{Txc,PsPot_UPF},
    psi::AbstractArray{ComplexF64},
    Hpsi::AbstractArray{ComplexF64}
) where Txc <: AbstractXCCalculator
    
    ik = Ham.ik
    Nstates = size(psi, 2)

    Nspecies = Ham.atoms.Nspecies
    Natoms = Ham.atoms.Natoms
    atm2species = Ham.atoms.atm2species

    NbetaNL = Ham.pspotNL.NbetaNL
    Deeq = Ham.pspotNL.Deeq
    nh = Ham.pspotNL.nh
    indv_ijkb0 = Ham.pspotNL.indv_ijkb0
    betaNL_k = Ham.pspotNL.betaNL[ik]

    # Use * directly?
    betaNL_psi = calc_betaNL_psi(ik, Ham.pspotNL, psi)

    ps = zeros(ComplexF64, NbetaNL, Nstates)

    for isp in 1:Nspecies

        if nh[isp] == 0
            continue
        end

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