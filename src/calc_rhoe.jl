function calc_rhoe!( Ham::Hamiltonian, psiks::BlochWavefunc, Rhoe::Array{Float64,2}; renormalize=true)

    pw = Ham.pw
    Focc = Ham.electrons.Focc
    Nspin = Ham.electrons.Nspin
    Nelectrons_true = Ham.electrons.Nelectrons

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
    for i in 1:length(Rhoe)
        if Rhoe[i] < eps()
            Rhoe[i] = eps()
        end
    end

    # renormalize
    if renormalize
        integ_rho = sum(Rhoe)*CellVolume/Npoints
        for ip in eachindex(Rhoe)
            Rhoe[ip] = Nelectrons_true/integ_rho * Rhoe[ip]
        end
    end

    # Symmetrize Rhoe if needed
    if Ham.sym_info.Nsyms > 1
        symmetrize_rhoe!( Ham.pw, Ham.sym_info, Ham.rhoe_symmetrizer, Rhoe )
    end

    return
end

function calc_rhoe( Ham::Hamiltonian, psiks::BlochWavefunc )
    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Rhoe = zeros(Float64, Npoints, Nspin)
    calc_rhoe!( Ham, psiks, Rhoe )
    return Rhoe
end


