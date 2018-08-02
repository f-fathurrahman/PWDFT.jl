function calc_rhoe( pw::PWGrid, Focc::Array{Float64,2},
                    psik::Array{Array{ComplexF64,2},1}; renormalize=true )
    CellVolume  = pw.CellVolume
    Ns = pw.Ns
    Nkpt = pw.gvecw.kpoints.Nkpt
    Ngw = pw.gvecw.Ngw
    wk = pw.gvecw.kpoints.wk
    Npoints = prod(Ns)
    Nstates = size(psik[1])[2]

    cpsi = zeros(ComplexF64, Npoints, Nstates)
    psiR = zeros(ComplexF64, Npoints, Nstates)
    rho = zeros(Float64,Npoints)

    for ik = 1:Nkpt
        cpsi[:,:] .= 0.0 + im*0.0
        # Transform to real space
        idx = pw.gvecw.idx_gw2r[ik]
        psi = psik[ik]
        cpsi[idx,:] = psi[:,:]
        psiR = G_to_R(pw, cpsi)
        # orthonormalization in real space
        ortho_gram_schmidt!( Nstates, psiR )
        psiR = sqrt(Npoints/CellVolume)*psiR
        #
        for ist = 1:Nstates
            for ip = 1:Npoints
                rho[ip] = rho[ip] + wk[ik]*Focc[ist,ik]*real( conj(psiR[ip,ist])*psiR[ip,ist] )
            end
        end
    end

    # Ensure that there is no negative rhoe
    for ip = 1:Npoints
        if rho[ip] < eps()
            rho[ip] = eps()
        end
    end

    # renormalize
    if renormalize
        integ_rho = sum(rho)*CellVolume/Npoints
        Nelectrons = sum(Focc)/Nkpt
        rho = Nelectrons/integ_rho * rho
    end

    return rho
end


"""
Calculate electron density from one kpoint.
"""
function calc_rhoe( ik::Int64, pw::PWGrid, Focc::Array{Float64,2},
                    psi::Array{ComplexF64,2}; renormalize=true )
    CellVolume  = pw.CellVolume
    Ns = pw.Ns
    Npoints = prod(Ns)
    Nstates = size(psi)[2]

    # Transform to real space
    cpsi = zeros( ComplexF64, Npoints, Nstates )
    idx = pw.gvecw.idx_gw2r[ik]
    cpsi[idx,:] = psi[:,:]
    psiR = G_to_R(pw, cpsi)

    # orthonormalization in real space
    ortho_gram_schmidt!( Nstates, psiR )
    psiR = sqrt(Npoints/CellVolume)*psiR

    rho = zeros(Float64,Npoints)
    for ist = 1:Nstates
        for ip = 1:Npoints
            rho[ip] = rho[ip] + Focc[ist,ik]*real( conj(psiR[ip,ist])*psiR[ip,ist] )
        end
    end

    # Ensure that there is no negative rhoe
    for ip = 1:Npoints
        if rho[ip] < eps()
            rho[ip] = eps()
        end
    end
    # renormalize
    if renormalize
        integ_rho = sum(rho)*CellVolume/Npoints
        Nelectrons = sum(Focc[:,ik])
        rho = Nelectrons/integ_rho * rho
    end

    return rho
end



# special case of non-spin pol and only 1 kpoint
function calc_rhoe( Ham::Hamiltonian,
                    psi::Array{ComplexF64,2}; renormalize=true )
    
    @assert( Ham.ik == 1 )
    @assert( Ham.ispin == 1 )

    ik = 1
    ispin = 1

    pw = Ham.pw
    Focc = Ham.electrons.Focc

    CellVolume  = pw.CellVolume
    Ns = pw.Ns
    Npoints = prod(Ns)
    Nstates = size(psi)[2]

    # Transform to real space
    cpsi = zeros( ComplexF64, Npoints, Nstates )
    idx = pw.gvecw.idx_gw2r[ik]
    cpsi[idx,:] = psi[:,:]
    psiR = G_to_R(pw, cpsi)

    # orthonormalization in real space
    ortho_gram_schmidt!( Nstates, psiR )
    psiR = sqrt(Npoints/CellVolume)*psiR

    Rhoe = zeros(Float64,Npoints)
    for ist = 1:Nstates
        for ip = 1:Npoints
            Rhoe[ip] = Rhoe[ip] + Focc[ist,ik]*real( conj(psiR[ip,ist])*psiR[ip,ist] )
        end
    end

    # Ensure that there is no negative rhoe
    for rho in Rhoe
        if rho < eps()
            rho = eps()
        end
    end
    
    # renormalize
    if renormalize
        integ_Rhoe = sum(Rhoe)*CellVolume/Npoints
        Nelectrons = sum(Focc[:,ik])
        Rhoe = Nelectrons/integ_Rhoe * Rhoe
    end

    return Rhoe
end