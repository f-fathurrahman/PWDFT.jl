function PAW_atomic_becsum!( Ham::Hamiltonian; starting_magnetization=nothing )
    PAW_atomic_becsum!(Ham.atoms, Ham.pspots, Ham.pspotNL,
        Nspin=Ham.electrons.Nspin, starting_magnetization=starting_magnetization)
end


# For compatibility
function PAW_atomic_becsum(atoms, pspots, pspotNL; Nspin=1, starting_magnetization=nothing)
    PAW_atomic_becsum!(atoms, pspots, pspotNL,
        Nspin=Nspin, starting_magnetization=starting_magnetization)
    return copy(pspotNL.becsum) # the copy is returned as it may become different from pspotNL.becsum
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
    starting_magnetization=nothing,
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

    @info "PAW_atomic_becsum: Nspin=$(Nspin)"
    # No starting magnetization is give, set them to a default value
    if (Nspin == 2) && isnothing(starting_magnetization)
        starting_magnetization = 0.1*ones(Nspecies)
        @info "Using default starting magnetization = $(starting_magnetization)"
    end

    becsum = pspotNL.paw.becsum # use PAW becsum
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

    # Need to symmetryize?


    println("becsum[1,1,1] after = ", becsum[1,1,1])
    println("becsum[2,1,1] after = ", becsum[2,1,1])
    if Nspin == 2
        println("becsum[1,1,2] after = ", becsum[1,1,2])
        println("becsum[2,1,2] after = ", becsum[2,1,2])
    end

    return
end
