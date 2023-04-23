# Sum up pfuncs x occupation to build radial density's angular momentum components.
# Based on PAW_rho_lm in paw_onecenter.f90 of QE.
#
function PAW_rho_lm!(
    AE::Bool, # if not AE it is PS
    ia::Int64,
    atoms::Atoms, pspots, pspotNL,
    becsum, rho_lm
)

    atm2species = atoms.atm2species
    isp = atm2species[ia]

    if AE
        pfunc = pspots[isp].paw_data.pfunc
    else
        pfunc = pspots[isp].paw_data.ptfunc
    end
    aug = pspots[isp].qfuncl # only used for is AE == false
  
    Nrmesh = pspots[isp].Nr
    nh = pspotNL.nh
    indv = pspotNL.indv
    lpx = pspotNL.lpx
    nhtolm = pspotNL.nhtolm
    lpl = pspotNL.lpl
    ap = pspotNL.ap

    Nspin = size(becsum,3)

    # initialize density
    fill!(rho_lm, 0.0)
    
    for ispin in 1:Nspin
        ijh = 0 
        # loop on all pfunc for this kind of pseudo 
        for ih in 1:nh[isp], jh in ih:nh[isp]
            ijh = ijh + 1 
            nb = indv[ih,isp] 
            mb = indv[jh,isp] 
            nmb = floor(Int64, (mb*(mb-1))/2 + nb)  # mb has to be >= nb 
            if abs(becsum[ijh,ia,ispin]) < 1e-12
                continue
            end 
            # Loop over angular momentum angular_momentum
            for lp in 1:lpx[nhtolm[jh,isp], nhtolm[ih,isp]] # lmaxq**2 
                # the lpl array contains the possible combination of LM,lm_j,lm_j that 
                # have non-zero a_{LM}^{(lm)_i(lm)_j} (it saves some loops) 
                lm = lpl[nhtolm[jh,isp], nhtolm[ih,isp], lp] 
                # 
                # becsum already contains a factor 2 for off-diagonal pfuncs 
                pref = becsum[ijh,ia,ispin] * ap[lm, nhtolm[ih,isp], nhtolm[jh,isp]] 
                #
                rho_lm[1:Nrmesh,lm,ispin] .+= pref * pfunc[1:Nrmesh,nb,mb] 
                if !AE
                    # if I'm doing the pseudo part I have to add the augmentation charge 
                    l = floor(Int64, sqrt(lm-1) ) # l has to start from zero, lm = l**2 + m
                    rho_lm[1:Nrmesh,lm,ispin] .+= pref * aug[1:Nrmesh,nmb,l+1]
                    # XXX: qfuncl is originally indexed from 0, we offset it by 1 here.
                end
            end # lp  
        end # ih, jh
    end
    return
end