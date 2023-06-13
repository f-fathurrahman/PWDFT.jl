function PAW_atomic_becsum!( Ham::Hamiltonian )
    PAW_atomic_becsum!(Ham.atoms, Ham.pspots, Ham.pspotNL, Nspin=Ham.electrons.Nspin)
end


# For compatibility
function PAW_atomic_becsum(atoms, pspots, pspotNL; Nspin=1)
    PAW_atomic_becsum!(atoms, pspots, pspotNL, Nspin=Nspin)
    return copy(pspotNL.becsum) # the copy is returned
end


# Initialize becsum with atomic occupations (for PAW atoms only).  
# NOTE: requires exact correspondence chi <--> beta in the atom,
# that is that all wavefunctions considered for PAW generation are
# counted in chi (otherwise the array "oc" does not correspond to beta).
#
# Note that symmetrization of becsum is not handled in this function.
# pspotNL.becsum is modified
function PAW_atomic_becsum!(
    atoms::Atoms,
    pspots::Vector{PsPot_UPF},
    pspotNL::PsPotNL_UPF;
    Nspin=1
)

    # Early return if no PAW species are present
    ok_paw = any(pspotNL.are_paw)
    if !ok_paw
        return
    end

    nhm = pspotNL.nhm
    nh = pspotNL.nh
    indv = pspotNL.indv
    nhtol = pspotNL.nhtol

    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species

    starting_magnetization = zeros(Nspecies)

    becsum = pspotNL.becsum
    fill!(becsum, 0.0)
  
    #!
    #! Add a bit of random noise if not starting from atomic or saved wfcs:
    #IF ( starting_wfc=='atomic+random') noise = 0.05_DP
    #IF ( starting_wfc=='random')        noise = 0.10_DP
    #noise = 0.1
    noise = 0.0
      
    for ia in 1:Natoms
        isp = atm2species[ia]
        psp = pspots[isp]
        if Nspin == 2
            μ = starting_magnetization[isp]
        end
        if psp.is_paw
            ijh = 1
            for ih in 1:nh[isp]
                nb = indv[ih,isp]
                ocnb = psp.paw_data.oc[nb]
                #
                if Nspin == 1
                    becsum[ijh,ia,1] = ocnb / (2*nhtol[ih,isp] + 1)
                #
                else # if Nspin == 2
                    becsum[ijh,ia,1] = 0.5*(1.0 + μ)*ocnb / (2*nhtol[ih,isp] + 1)
                    becsum[ijh,ia,2] = 0.5*(1.0 - μ)*ocnb / (2*nhtol[ih,isp] + 1)
                end
                #
                ijh = ijh + 1
                for jh in (ih + 1):nh[isp]
                    #mb = indv(jh,nt)
                    for ispin in 1:Nspin
                        if noise > 0.0
                            becsum[ijh,ia,ispin] += noise * 2.0 * ( 0.5 - rand() )
                        end
                    end
                    ijh = ijh + 1
                end # jh loop
            end # ih loop
        end # if psp.is_paw
    end # ia loop

    return
end
