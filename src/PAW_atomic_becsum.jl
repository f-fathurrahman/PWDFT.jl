
# Initialize becsum with atomic occupations (for PAW atoms only).  
# NOTE: requires exact correspondence chi <--> beta in the atom,
# that is that all wavefunctions considered for PAW generation are
# counted in chi (otherwise the array "oc" does not correspond to beta).
function PAW_atomic_becsum(atoms, pspots, pspotNL; Nspin=1)

    # 
    # REAL(DP), INTENT(INOUT) :: becsum(nhm*(nhm+1)/2,nat,nspin)

    nhm = pspotNL.nhm
    nh = pspotNL.nh
    indv = pspotNL.indv
    nhtol = pspotNL.nhtol

    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species

    starting_magnetization = zeros(Nspecies)

    Nbecsum = Int64( nhm * (nhm + 1)/2 )
    becsum = zeros(Float64, Nbecsum, Natoms, Nspin)


    # IF (.NOT. okpaw) RETURN
  
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
                elseif Nspin == 2
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

    println("sum becsum = ", sum(becsum))

    #!
    #! ... copy becsum in scf structure and symmetrize it
    #rho%bec(:,:,:) = becsum(:,:,:) 
    #CALL PAW_symmetrize( rho%bec )

    return becsum
end
