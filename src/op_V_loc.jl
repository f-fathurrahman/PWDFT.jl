function op_V_loc( pw::PWGrid, V_loc::Array{Float64,1}, psi::Array{Complex128,2} )
    Ns = pw.Ns
    Ω  = pw.Ω
    Npoints = prod(Ns)
    Nstates = size(psi)[2]

    ctmp = zeros(Complex128, Npoints, Nstates)
    idx = pw.gvecw.idx_gw2r
    for ic = 1:Nstates
        ctmp[idx,ic] = psi[:,ic]
    end

    # get values of psi in real space grid via forward transform
    ctmp = G_to_R( pw, ctmp )

    for ist = 1:Nstates
        for ip = 1:Npoints
            ctmp[ip,ist] = V_loc[ip]*ctmp[ip,ist]
        end
    end

    cVpsi = R_to_G( pw, ctmp )
    return cVpsi[idx,:]
end

function op_V_loc( Ham::PWHamiltonian, psi::Array{Complex128,2} )
    V_loc = Ham.potentials.Ps_loc + Ham.potentials.Hartree + Ham.potentials.XC
    return op_V_loc( Ham.pw, V_loc, psi )
end

function op_V_Ps_loc( Ham::PWHamiltonian, psi::Array{Complex128,2} )
    return op_V_loc( Ham.pw, Ham.potentials.Ps_loc, psi )
end

#
# single-column version
#
function op_V_loc( pw::PWGrid, V_loc::Array{Float64,1}, psi::Array{Complex128,1} )
    Ns = pw.Ns
    Ω  = pw.Ω
    Npoints = prod(Ns)

    ctmp = zeros(Complex128, Npoints)
    idx = pw.gvecw.idx_gw2r
    ctmp[idx] = psi[:]

    # get values of psi in real space grid via forward transform
    ctmp = G_to_R( pw, ctmp )

    for ip = 1:Npoints
        ctmp[ip] = V_loc[ip]*ctmp[ip]
    end

    cVpsi = R_to_G( pw, ctmp )
    return cVpsi[idx]
end

function op_V_loc( Ham::PWHamiltonian, psi::Array{Complex128,1} )
    V_loc = Ham.potentials.Ps_loc + Ham.potentials.Hartree + Ham.potentials.XC
    return op_V_loc( Ham.pw, V_loc, psi )
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
    ctmp = G_to_R( pw, ctmp )

    cVpsi = R_to_G( pw, V_Ps_loc .* ctmp )
    return cVpsi[idx]
end
