function calc_rhoe(
    Nelectrons_true::Float64,
    pw::PWGrid, Focc::Array{Float64,2},
    psiks::BlochWavefunc, Nspin::Int64,

    ;
    renormalize=true,
)

    CellVolume  = pw.CellVolume
    Ns = pw.Ns
    Nkpt = pw.gvecw.kpoints.Nkpt
    Ngw = pw.gvecw.Ngw
    wk = pw.gvecw.kpoints.wk
    Npoints = prod(Ns)
    Nstates = size(psiks[1])[2]

    #cpsi = zeros(ComplexF64, Npoints, Nstates)
    psiR = zeros(ComplexF64, Npoints, Nstates)
    Rhoe = zeros(Float64, Npoints, Nspin)

    for ispin = 1:Nspin
    for ik = 1:Nkpt

        ikspin = ik + (ispin - 1)*Nkpt

        #cpsi[:,:] .= 0.0 + im*0.0
        psiR .= 0.0 + im*0.0
        
        # Transform to real space
        idx = pw.gvecw.idx_gw2r[ik]
        psi = psiks[ikspin]
        
        # Using the following loop gives not too much improvement
        # over slicing-operations
        #for ist = 1:Nstates
        #    for ig = 1:Ngw[ik]
        #        psiR[idx[ig],ist] = psi[ig,ist]
        #    end
        #end

        psiR[idx,:] = psi[:,:]
        G_to_R!(pw, psiR)

        # orthonormalization in real space
        ortho_sqrt!( psiR )
        psiR = sqrt(Npoints/CellVolume)*psiR
        #
        for ist = 1:Nstates
            w = wk[ik]*Focc[ist,ikspin]
            for ip = 1:Npoints
                Rhoe[ip,ispin] = Rhoe[ip,ispin] + w*real( conj(psiR[ip,ist])*psiR[ip,ist] )
            end
        end
    end # ik
    end # ikspin

    # Ensure that there is no negative rhoe
    for rho in Rhoe
        if rho < eps()
            rho = 0.0
        end
    end

    # renormalize
    if renormalize
        integ_rho = sum(Rhoe)*CellVolume/Npoints
        Rhoe = Nelectrons_true/integ_rho * Rhoe
    end

    return Rhoe
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
    ortho_sqrt!( psiR )
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
    ortho_sqrt!( Nstates, psiR )
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
