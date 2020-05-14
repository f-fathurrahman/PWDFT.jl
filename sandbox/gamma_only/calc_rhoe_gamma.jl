import PWDFT: calc_rhoe!, calc_rhoe

function calc_rhoe!(
    Ham::HamiltonianGamma, psis::BlochWavefuncGamma, Rhoe::Array{Float64,2}; renormalize=false
)

    pw = Ham.pw
    Focc = Ham.electrons.Focc
    Nspin = Ham.electrons.Nspin
    Nelectrons_true = Ham.electrons.Nelectrons

    CellVolume  = pw.CellVolume
    Ns = pw.Ns
    Ngw = pw.gvecw.Ngw
    Npoints = prod(Ns)
    Nstates = size(psis.data[1],2)

    psiR = zeros(ComplexF64, Npoints, Nstates)

    # dont forget to zero out the Rhoe first
    Rhoe[:,:] .= 0.0

    idx_gw2r = pw.gvecw.idx_gw2r
    idx_gw2rm = pw.gvecw.idx_gw2rm

    for ispin = 1:Nspin

        psiR .= 0.0 + im*0.0
        
        # Transform to real space
        for ist in 1:Nstates
            psiR[1,ist] = psis.data[ispin][1,ist]
            for igw in 2:Ngw
                ip = idx_gw2r[igw]
                ipm = idx_gw2rm[igw]
                psiR[ip,ist] = psis.data[ispin][igw,ist]
                psiR[ipm,ist] = conj(psis.data[ispin][igw,ist])
            end
        end
        
        G_to_R!(pw, psiR)

        # orthonormalization in real space
        #ortho_sqrt!( psiR )
        ##ortho_gram_schmidt!( psiR ) # for comparison
        #psiR = sqrt(Npoints/CellVolume)*psiR
        
        # by pass orthonormalization, only use scaling
        psiR = sqrt(Npoints/CellVolume)*sqrt(Npoints)*psiR

        for ist = 1:Nstates, ip = 1:Npoints
            Rhoe[ip,ispin] = Rhoe[ip,ispin] + Focc[ist,ispin]*real( conj(psiR[ip,ist])*psiR[ip,ist] )
        end
    end # ispin

    # Ensure that there is no negative rhoe
    #for i in 1:length(Rhoe)
    #    if Rhoe[i] < eps()
    #        Rhoe[i] = eps()
    #    end
    #end

    # renormalize
    if renormalize
        integ_rho = sum(Rhoe)*CellVolume/Npoints
        for ip in eachindex(Rhoe)
            Rhoe[ip] = Nelectrons_true/integ_rho * Rhoe[ip]
        end
    end

    return
end

# TODO: Remove type annotation
function calc_rhoe( Ham::HamiltonianGamma, psis::BlochWavefuncGamma; renormalize=false )
    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Rhoe = zeros(Float64, Npoints, Nspin)
    calc_rhoe!( Ham, psis, Rhoe, renormalize=renormalize )
    return Rhoe
end


