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

    rho = zeros(Float64,Npoints)
    for ist = 1:Nstates
        for ip = 1:Npoints
            rho[ip] = rho[ip] + Focc[ist]*real( conj(psiR[ip,ist])*psiR[ip,ist] )
        end
    end

    # Ensure that there is no negative rhoe
    for ip = 1:Nstates
        if rho[ip] < eps()
            rho[ip] = eps()
        end
    end
    # renormalize
    integ_rho = sum(rho)*Ω/Npoints
    Nelectrons = sum(Focc)
    rho = Nelectrons/integ_rho * rho

    return rho
end
