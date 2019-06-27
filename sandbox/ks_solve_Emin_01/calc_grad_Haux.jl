function calc_grad_Haux( Ham::Hamiltonian, psi::Array{ComplexF64,2}, eta::Array{Float64,1}, kT::Float64 )

    ik = Ham.ik
    ispin = Ham.ispin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    ikspin = ik + (ispin - 1)*Nkpt

    f = @view Ham.electrons.Focc[:,ikspin]

    #println("f = ", f)
    
    Ngw     = size(psi)[1]
    Nstates = size(psi)[2]

    # gradients
    g_psi = zeros(ComplexF64, Ngw, Nstates)
    g_Haux = zeros(ComplexF64, Nstates, Nstates)

    # primary search direction
    Δ_psi = zeros(ComplexF64, Ngw, Nstates)
    Δ_Haux = zeros(ComplexF64, Nstates, Nstates)

    Hpsi = op_H( Ham, psi )
    Hsub = psi' * Hpsi

    for ist = 1:Nstates
        Δ_psi[:,ist] = Hpsi[:,ist]
        for jst = 1:Nstates
            Δ_psi[:,ist] = Δ_psi[:,ist] - Hsub[jst,ist]*psi[:,jst]
        end
        g_psi[:,ist] = f[ist]*Δ_psi[:,ist]
    end

    df_deta = zeros(Nstates)
    for ist = 1:Nstates
        df_deta[ist] = -0.5*f[ist]*(1.0 - 0.5*f[ist])/kT  # factor 0.5 for non-spinpol
    end

    ss = 0.0
    for ist = 1:Nstates
        ss = ss + 0.5*f[ist]*(1.0 - 0.5*f[ist])
    end
    dmu_deta = zeros(Nstates)
    for ist = 1:Nstates
        dmu_deta[ist] = 0.5*f[ist]*(1.0 - 0.5*f[ist])/ss
    end

    #for ist = 1:Nstates
    #    @printf("dmu_deta = %18.10f, df_deta= %18.10f\n", dmu_deta[ist], df_deta[ist])
    #end


    Δ_Haux[:,:] = copy(Hsub)

    # diagonal
    for ist = 1:Nstates
        term1 = ( Hsub[ist,ist] - eta[ist] ) * df_deta[ist]
        term2 = 0.0 + im*0.0
        for kst = 1:Nstates
            term2 = term2 + (Hsub[kst,kst] - eta[kst])*df_deta[kst]
        end
        g_Haux[ist,ist] = term1 - dmu_deta[ist]*term2
        Δ_Haux[ist,ist] = Hsub[ist,ist] - eta[ist]
    end

    # off diagonal
    for ist = 1:Nstates
        for jst = (ist+1):Nstates
            g_Haux[ist,jst] = Hsub[ist,jst] * (f[ist] - f[jst]) / (eta[ist] - eta[jst])
            g_Haux[jst,ist] = g_Haux[ist,jst]
        end
    end

    return g_psi, g_Haux, Δ_psi, Δ_Haux
end
