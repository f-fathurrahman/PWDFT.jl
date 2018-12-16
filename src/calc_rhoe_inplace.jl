function calc_rhoe!( Ham::Hamiltonian, psiks::BlochWavefunc, Rhoe::Array{Float64,2} )
    calc_rhoe!(
        Rhoe,
        Ham.electrons.Nelectrons,
        Ham.pw, Ham.electrons.Focc, 
        psiks, Ham.electrons.Nspin )
end

# in-place version
function calc_rhoe!(
    Rhoe::Array{Float64,2},
    Nelectrons_true::Float64,
    pw::PWGrid, Focc::Array{Float64,2},
    psiks::BlochWavefunc, Nspin::Int64;
    renormalize=true
)

    CellVolume  = pw.CellVolume
    Ns = pw.Ns
    Nkpt = pw.gvecw.kpoints.Nkpt
    Ngw = pw.gvecw.Ngw
    wk = pw.gvecw.kpoints.wk
    Npoints = prod(Ns)
    Nstates = size(psiks[1])[2]

    psiR = zeros(ComplexF64, Npoints, Nstates)

    # dont forget to zero out the Rhoe first
    Rhoe[:,:] .= 0.0

    for ispin = 1:Nspin
    for ik = 1:Nkpt

        ikspin = ik + (ispin - 1)*Nkpt

        #cpsi[:,:] .= 0.0 + im*0.0
        psiR .= 0.0 + im*0.0
        
        # Transform to real space
        idx = pw.gvecw.idx_gw2r[ik]
        psi = psiks[ikspin]

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
            rho = eps()
        end
    end

    # renormalize
    if renormalize
        integ_rho = sum(Rhoe)*CellVolume/Npoints
        Rhoe = Nelectrons_true/integ_rho * Rhoe
    end

    return
end