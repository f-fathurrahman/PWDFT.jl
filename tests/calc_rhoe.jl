function calc_rhoe( pw::PWGrid, Focc, psi::Array{Complex128,2} )
    Ω  = pw.Ω
    Ns = pw.Ns
    Npoints = prod(Ns)
    Nstates = size(psi)[2]

    # Transform to real space
    cpsi = zeros( Complex128, Npoints, Nstates )
    cpsi[pw.gvecw.idx_gw2r,:] = psi[:,:]
    psiR = G_to_R(Ns, cpsi)

    # orthonormalization in real space
    ortho_gram_schmidt!( Nstates, psiR )
    scale!( sqrt(Npoints/Ω), psiR )

    rho = Array{Float64}(Npoints)
    for is = 1:Nstates
        for ip = 1:Npoints
            rho[ip] = rho[ip] + Focc[is]*real( conj(psiR[ip,is])*psiR[ip,is] )
        end
    end

    # Ensure that there is no negative rhoe
    for ip = 1:Nstates
        if rho[ip] < 0.0
            rho[ip] = eps()
        end
    end
    # renormalize
    integ_rho = sum(rho)*Ω/Npoints
    Nelectrons = sum(Focc)
    rho = Nelectrons/integ_rho * rho

    return rho
end
