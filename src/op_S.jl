function op_S(
    Ham::Hamiltonian{Txc,PsPot_UPF},
    psi::AbstractArray{ComplexF64}
) where Txc <: AbstractXCCalculator

    Spsi = zeros(ComplexF64, size(psi))
    op_S!(Ham, psi, Spsi)
    return Spsi
end


function op_S!(
    Ham::Hamiltonian{Txc,PsPot_UPF},
    psi::AbstractArray{ComplexF64},
    Spsi::AbstractArray{ComplexF64}
) where Txc <: AbstractXCCalculator

    # Check Vnl_KB construction
    ik = Ham.ik
    ispin = Ham.ispin # XXX: used?

    pspotNL = Ham.pspotNL
    pspots = Ham.pspots
    atoms = Ham.atoms

    Ngwk = size(psi,1) # check againsts Ngw[k]
    Nstates = size(psi,2)

    Vnl_KB = pspotNL.betaNL[ik]

    betaNL_psi = Vnl_KB' * psi  # XXX precalculate this

    Nspecies = atoms.Nspecies
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    nkb = pspotNL.NbetaNL
    nh = pspotNL.nh
    indv_ijkb0 = pspotNL.indv_ijkb0
    qq_at = pspotNL.qq_at

    ps = zeros(ComplexF64, nkb, Nstates)
    
    @views Spsi[:,:] = psi[:,:]

    for isp in 1:Nspecies
        #
        if pspots[isp].is_ultrasoft
            for ia in 1:Natoms
                if atm2species[ia] != isp
                    continue
                end
                idx1 = (indv_ijkb0[ia]+1):(indv_ijkb0[ia]+nh[isp])
                ps[idx1,1:Nstates] = qq_at[1:nh[isp],1:nh[isp],ia] * betaNL_psi[idx1,:]

            end
        else
            if nh[isp] > 0
                for ia in 1:Natoms
                    if atm2species[ia] != isp
                        continue
                    end
                    idx1 = (indv_ijkb0[ia]+1):(indv_ijkb0[ia]+nh[isp])
                    ps[idx1,1:Nstates] .= 0.0 + im*0.0
                end
            end
        end
    end

    @views Spsi[:,:] += Vnl_KB * ps

    return

end