function dot_ElecGradient( v1::ElecGradient, v2::ElecGradient )
    Nkspin = length(v1.psiks)
    ss = 0.0
    for i in 1:Nkspin
        ss = ss + 2.0*real( dot(v1.psiks[i], v2.psiks[i]) )
        ss = ss + 2.0*real( dot(v1.Haux[i], v2.Haux[i]) )
    end
    return ss
end

function compute!(
    Ham::Hamiltonian,
    evars::ElecVars,
    g::ElecGradient, Kg::ElecGradient, kT::Float64,
    rotPrevCinv, rotPrev
)

    Etot = calc_energies_grad!( Ham, evars, g, Kg, kT )

    Nkspin = length(evars.psiks)

    for i in 1:Nkspin
        g.psiks[i] = g.psiks[i] * rotPrevCinv[i]
        Kg.psiks[i] = Kg.psiks[i] * rotPrevCinv[i]
        g.Haux[i] = rotPrev[i] * g.Haux[i] * rotPrev[i]'
        Kg.Haux[i] = rotPrev[i] * Kg.Haux[i] * rotPrev[i]'
    end

    # No caching is done (for SubspaceRotationAdjutst)

    return Etot

end

function do_step!(
    α::Float64, evars::ElecVars, d::ElecGradient,
    rotPrev, rotPrevC, rotPrevCinv
)
    
    Nkspin = length(psiks)
    Nstates = size(psiks[1],2)
    
    Haux = zeros(ComplexF64,Nstates,Nstates)
    rot = zeros(ComplexF64,Nstates,Nstates)
    rotC = zeros(ComplexF64,Nstates,Nstates)

    for i in 1:Nkspin
        evars.psiks[i] = evars.psiks[i] + α*d.psiks[i]

        # Haux fillings:
        Haux = diagm( 0 => eVars.Haux_eigs[:,i] )
        
        #axpy(alpha, rotExists ? dagger(rotPrev[q])*dir.Haux[q]*rotPrev[q] : dir.Haux[q], Haux);
        Haux = Haux + α*( rotPrev[i]' * d.Haux[i] * rotPrev[i] )
        
        #Haux.diagonalize(rot, eVars.Haux_eigs[q]); //rotation chosen to diagonalize auxiliary matrix
        evars.Haux_eigs[:,i], rot = eigen(Haux)

        #rotC = rot
        #eVars.orthonormalize(q, &rotC);
        Udagger = inv( sqrt( evars.psiks[i]' * evars.psiks[i] ) )
        rotC = Udagger*rot
        evars.psiks[i] = evars.psiks[i]*rotC
        
        rotPrev[i] = rotPrev[i] * rot
        rotPrevC[i] = rotPrevC[i] * rotC
        rotPrevCinv[i] = inv(rotC) * rotPrevCinv[q]

    end
    
    return 
end

function constrain_search_dir!( d::ElecGradient, evars::ElecVars )
    Nkspin = length(evars.psiks)
    for i in 1:Nkspin
        d.psiks[i] = d.psiks[i] - evars.psiks[i] * ( evars.psiks[i]' * d.psiks[i] )
    end
    return
end
