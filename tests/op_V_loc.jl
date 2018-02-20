function op_V_Ps_loc( Ham::PWHamiltonian, psi::Array{Complex128,2} )
    #
    pw = Ham.pw
    #
    Ns = pw.Ns
    Ω  = pw.Ω
    Npoints = prod(Ns)
    Nstates = size(psi)[2]
    #
    V_Ps_loc = Ham.potentials.Ps_loc

    ctmp = zeros(Complex128, Npoints, Nstates)
    idx = pw.gvecw.idx_gw2r
    for ic = 1:Nstates
        ctmp[idx,ic] = psi[:,ic]
    end

    # get values of psi in real space grid via forward transform
    ctmp = G_to_R( Ns, ctmp )

    for ist = 1:Nstates
        for ip = 1:Npoints
            ctmp[ip,ist] = V_Ps_loc[ip]*ctmp[ip,ist]
        end
    end

    cVpsi = R_to_G( Ns, ctmp )
    return cVpsi[idx,:]
end

function op_V_Ps_loc( Ham::PWHamiltonian, psi::Array{Complex128,1} )
    #
    pw = Ham.pw
    #
    Ns = pw.Ns
    Ω  = pw.Ω
    Npoints = prod(Ns)
    #
    V_Ps_loc = Ham.potentials.Ps_loc
    #
    ctmp = zeros(Complex128, Npoints)
    idx = pw.gvecw.idx_gw2r
    ctmp[idx] = psi[:]
    #
    # get values of psi in real space grid via forward transform
    ctmp = G_to_R( Ns, ctmp )

    cVpsi = R_to_G( Ns, V_Ps_loc .* ctmp )
    return cVpsi[idx]
end
