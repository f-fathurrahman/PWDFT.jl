function update_occ!( Ham, evars::ElecVars, kT )
    
    Nspin = Ham.electrons.Nspin
    Nelectrons = Ham.electrons.Nelectrons
    wk = Ham.pw.gvecw.kpoints.wk

    Ham.electrons.ebands = copy(evars.Haux_eigs)

    Focc, E_fermi = calc_Focc( Nelectrons, wk, kT, Ham.electrons.ebands, Nspin )
    mTS = calc_entropy( wk, kT, Ham.electrons.ebands, E_fermi, Nspin )
    
    Ham.electrons.Focc = copy(Focc)

    println("\nIn update_occ!")
    print_ebands(Ham.electrons, Ham.pw.gvecw.kpoints)


    return E_fermi, mTS

end


function calc_energies_grad!(
    Ham::Hamiltonian,
    evars::ElecVars,
    g::ElecGradient, Kg::ElecGradient, kT::Float64
)

    E_fermi, mTS = update_occ!( Ham, evars, kT )
    println("E_fermi = ", E_fermi)

    Rhoe = calc_rhoe( Ham, evars.psiks )
    update!( Ham, Rhoe )
    
    Ham.energies = calc_energies( Ham, evars.psiks )
    Ham.energies.mTS = mTS

    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nstates = Ham.electrons.Nstates
    Nkspin = Nkpt*Nspin

    #
    # Gradient for psiks
    #
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        i = ik + (ispin-1)*Nkpt
        calc_grad!( Ham, evars.psiks[i], g.psiks[i], evars.Hsub[i] )
        Kprec!( ik, Ham.pw, g.psiks[i], Kg.psiks[i] )
    end

    #
    # Gradient for Haux
    #
    fprime = zeros(Float64,Nstates)
    fprimeNum = zeros(Float64,Nstates)
    dmuNum = zeros(Float64,Nspin)
    dmuDen = zeros(Float64,Nspin)

    w = copy(Ham.pw.gvecw.kpoints.wk)
    if Nspin == 1
        w = w*2.0
    end

    for ispin in 1:Nspin, ik in 1:Nkpt
        i = ik + (ispin - 1)*Nkpt
        for ist in 1:Nstates
            fprime[ist] = smear_fermi_prime( evars.Haux_eigs[ist,i], E_fermi, kT )
            fprimeNum[ist] = fprime[ist] * ( real(evars.Hsub[i][ist,ist]) - evars.Haux_eigs[ist,i] )
        end
        dmuNum[ispin] = dmuNum[ispin] + w[ik] * sum(fprimeNum)
        dmuDen[ispin] = dmuDen[ispin] + w[ik] * sum(fprime)
    end

    dmuContrib = sum(dmuNum)/sum(dmuDen)
    dBzContrib = 0.0 # not used

    gradF0 = zeros(ComplexF64,Nstates,Nstates)
    gradF = zeros(ComplexF64,Nstates,Nstates)

    g_tmp = zeros(ComplexF64,Nstates,Nstates)
    for ispin in 1:Nspin, ik in 1:Nkpt
        #
        i = ik + (ispin - 1)*Nkpt
        #
        gradF0[:] = evars.Hsub[i] - diagm( 0 => evars.Haux_eigs[:,i] )
        gradF[:] = copy(gradF0)
        for ist in 1:Nstates
            gradF[ist,ist] = gradF0[ist,ist] - dmuContrib # FIXME: not tested for spinpol
        end
        g_tmp[:] = grad_smear( smear_fermi, smear_fermi_prime, evars.Haux_eigs[:,i], E_fermi, kT, gradF )
        g.Haux[i] = w[ik] * 0.5 * (g_tmp' + g_tmp)
        Kg.Haux[i] = -copy(gradF0) #-0.1*copy(gradF0)
    end

    return sum( Ham.energies )
end

function calc_grad!( Ham::Hamiltonian, psiks::BlochWavefunc, g::BlochWavefunc, Hsub )
    #
    Rhoe = calc_rhoe( Ham, psiks )
    update!( Ham, Rhoe )
    #
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    #
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        i = ik + (ispin-1)*Nkpt
        calc_grad!( Ham, psiks[i], g[i], Hsub[i] )
    end
    return
end

function calc_grad!( Ham::Hamiltonian, ψ::Array{ComplexF64,2}, g::Array{ComplexF64,2}, Hsub::Matrix{ComplexF64} )

    ik = Ham.ik
    ispin = Ham.ispin

    Nstates = size(ψ,2)
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    ikspin = ik + (ispin - 1)*Nkpt
    Focc = Ham.electrons.Focc
    wk_ik = Ham.pw.gvecw.kpoints.wk[ik]

    Hψ = op_H( Ham, ψ )
    Hsub[:] = ψ' * Hψ
    Hψ = Hψ - ψ*Hsub

    NkFull = prod(Ham.pw.gvecw.kpoints.mesh)
    for ist in 1:Nstates
        g[:,ist] = wk_ik * Focc[ist,ikspin] * Hψ[:,ist]
    end

    return

end

function Kprec!( ik::Int64, pw::PWGrid, psi::Array{ComplexF64,2}, Kpsi::Array{ComplexF64,2} )
    Ngw_ik  = size(psi)[1]
    Nstates = size(psi)[2]
    idx_gw2g = pw.gvecw.idx_gw2g[ik]
    G = pw.gvec.G
    k = pw.gvecw.kpoints.k[:,ik]
    for ist = 1:Nstates
        for igk = 1:Ngw_ik
            ig = idx_gw2g[igk]
            Gw2 = (G[1,ig] + k[1])^2 + (G[2,ig] + k[2])^2 + (G[3,ig] + k[3])^2
            Kpsi[igk,ist] = psi[igk,ist] / ( 1.0 + Gw2 )
        end
    end
    return
end

function dot_ElecGradient( v1::ElecGradient, v2::ElecGradient )
    Nkspin = length(v1.psiks)
    ss = 0.0
    for i in 1:Nkspin
        ss = ss + 2.0*real( dot(v1.psiks[i], v2.psiks[i]) )
        ss = ss + real( dot(v1.Haux[i], v2.Haux[i]) ) # no factor of 2
    end
    return ss
end

function dot_ElecGradient_v2( v1::ElecGradient, v2::ElecGradient )
    Nkspin = length(v1.psiks)
    ss = 0.0
    ss_Haux = 0.0
    for i in 1:Nkspin
        ss = ss + 2.0*real( dot(v1.psiks[i], v2.psiks[i]) )
        ss_Haux = ss_Haux + real( dot(v1.Haux[i], v2.Haux[i]) ) # no factor of 2
    end
    return ss, ss_Haux
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
    do_step!(α, α, evars, d, rotPrev, rotPrevC, rotPrevCinv)
    return
end

function do_step!(
    α::Float64, α_Haux::Float64, evars::ElecVars, d::ElecGradient,
    rotPrev, rotPrevC, rotPrevCinv
)
    
    Nkspin = length(evars.psiks)
    Nstates = size(evars.psiks[1],2)
    
    Haux = zeros(ComplexF64,Nstates,Nstates)
    rot = zeros(ComplexF64,Nstates,Nstates)
    rotC = zeros(ComplexF64,Nstates,Nstates)

    for i in 1:Nkspin
        evars.psiks[i] = evars.psiks[i] + α*d.psiks[i]*rotPrevC[i]

        # Haux fillings:
        Haux = diagm( 0 => evars.Haux_eigs[:,i] )
        
        #axpy(alpha, rotExists ? dagger(rotPrev[q])*dir.Haux[q]*rotPrev[q] : dir.Haux[q], Haux);
        Haux = Haux + α_Haux*( rotPrev[i]' * d.Haux[i] * rotPrev[i] )
        
        #Haux.diagonalize(rot, eVars.Haux_eigs[q]); //rotation chosen to diagonalize auxiliary matrix
        #evals, evecs = eigen(Haux)
        #println("evals = ", evals)
        evars.Haux_eigs[:,i], rot = eigen(Hermitian(Haux)) # need to symmetrize?
 
        #rotC = rot
        #eVars.orthonormalize(q, &rotC);
        Udagger = inv( sqrt( evars.psiks[i]' * evars.psiks[i] ) )
        rotC = Udagger*rot
        evars.psiks[i] = evars.psiks[i]*rotC
        
        rotPrev[i] = rotPrev[i] * rot
        rotPrevC[i] = rotPrevC[i] * rotC
        rotPrevCinv[i] = inv(rotC) * rotPrevCinv[i]

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
