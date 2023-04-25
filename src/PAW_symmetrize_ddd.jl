function PAW_symmetrize_ddd!( Ham::Hamiltonian )
    PAW_symmetrize_ddd!( Ham.atoms, Ham.pspots, Ham.sym_info, Ham.pspotNL,
        Ham.pspotNL.paw.ddd_paw )
    return
end

function PAW_symmetrize_ddd!(
    atoms::Atoms,
    pspots,
    sym_info::SymmetryInfo,
    pspotNL,
    ddd
)    
    #
    # => lm = l**2 + m
    # => ih = lm + (l+proj)**2  <-- if the projector index starts from zero!
    #       = lm + proj**2 + 2*l*proj
    #       = m + l**2 + proj**2 + 2*l*proj
    #        ^^^
    # Known ih and m_i I can compute the index oh of a different m = m_o but
    # the same augmentation channel (l_i = l_o, proj_i = proj_o):
    #  oh = ih - m_i + m_o
    # this expression should be general inside pwscf.

    Nsyms = sym_info.Nsyms
    if Nsyms == 1
        return
    end

    Nspin = size(ddd,3)
    
    D = Vector{Array{Float64,3}}(undef,4)
    D[1] = ones(1,1,Nsyms)
    D[2] = sym_info.D1
    D[3] = sym_info.D2
    D[4] = sym_info.D3

    atm2species = atoms.atm2species
    Natoms = atoms.Natoms

    nh = pspotNL.nh
    ijtoh = pspotNL.ijtoh
    nhtolm = pspotNL.nhtolm
    nhtol = pspotNL.nhtol

    irt = sym_info.irt
    dddsym = zeros(Float64, size(ddd))

    println("size ddd = ", size(ddd))

    usym = 1.0 / Nsyms
    
    for ispin in 1:Nspin, ia in 1:Natoms
        
        isp = atm2species[ia]

        # No need to symmetrize non-PAW atoms
        if !pspots[isp].is_paw
            continue
        end

        for ih in 1:nh[isp], jh in ih:nh[isp] # note: jh >= ih

            ijh = ijtoh[ih,jh,isp]
            
            lm_i = nhtolm[ih,isp]
            lm_j = nhtolm[jh,isp]
            
            l_i = nhtol[ih,isp]
            l_j = nhtol[jh,isp]
            
            m_i   = lm_i - l_i^2
            m_j   = lm_j - l_j^2
            
            for isym in 1:Nsyms
                ma = irt[isym,ia] 
                for m_o in 1:(2*l_i+1), m_u in 1:(2*l_j+1)
                    oh = ih - m_i + m_o
                    uh = jh - m_j + m_u
                    ouh = ijtoh[oh,uh,isp]
                    dddsym[ijh,ia,ispin] += D[l_i+1][m_o,m_i,isym]*D[l_j+1][m_u,m_j,isym]*usym*ddd[ouh,ma,ispin]
                end
            end
        end
    end

    #
    # Apply symmetrization:
    ddd[:,:,:] .= dddsym[:,:,:]

    return

end