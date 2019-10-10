"""
Prepare guess wavefunc
"""
function guess_evars!( Ham, evars; NiterMax=10 )

    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Rhoe = zeros(Float64,Npoints,Nspin)

    if Nspin == 1
        Rhoe[:,1] = guess_rhoe( Ham )
    else
        Rhoe = guess_rhoe_atomic( Ham, starting_magnetization=0.1*ones(Ham.atoms.Nspecies) )
    end
    update!(Ham, Rhoe)

    Ham.electrons.ebands[:,:] = diag_LOBPCG!( Ham, evars.ψ, verbose=true, NiterMax=NiterMax )
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        i = ik + (ispin - 1)*Nkpt
        evars.η[i] = diagm( 0 => Ham.electrons.ebands[:,i] )
    end

    return
end


"""
Check convergence in Haux
"""
function check_Hsub( Ham, evars )
    ψ = evars.ψ
    η = evars.η
    Nkspin = length(ψ)
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin
    Nstates = Ham.electrons.Nstates
    ss = 0.0
    for ispin in 1:Nspin, ik in 1:Nkpt
        #
        Ham.ik = ik
        Ham.ispin = ispin
        #
        i = ik + (ispin - 1)*Nkpt
        Hsub = Hermitian( ψ[i]' * (Ham*ψ[i]) )
        λ = eigvals(Hsub)
        for ist in 1:Nstates
            #@printf("Hsub-Haux: %d %18.10f\n", ist, abs(λ[ist] - real(η[i][ist,ist])))
            #ds = abs(real(Hsub[ist,ist] - η[i][ist,ist]))
            ds = abs(λ[ist] - real(η[i][ist,ist]))
            @printf("Hsub-Haux: %d %18.10f %18.10f %18.10f diff: %18.10f\n", ist, Ham.electrons.Focc[ist,i],
               λ[ist], real(η[i][ist,ist]), ds)
            ss = ss + ds
        end
    end
    println("ss = ", ss)
    return
end

"""
A convenience wrapper for `calc_Focc` and `calc_entropy`
Modify Ham.electrons.Focc and returns E_fermi and Entropy
"""
function set_occupations!( Ham, kT )
    Ham.electrons.Focc, E_fermi = calc_Focc(
        Ham.electrons.Nelectrons, Ham.pw.gvecw.kpoints.wk, kT, Ham.electrons.ebands, Ham.electrons.Nspin )
    Entropy = calc_entropy( Ham.pw.gvecw.kpoints.wk, kT,
        Ham.electrons.ebands, E_fermi, Ham.electrons.Nspin
    )
    return E_fermi, Entropy
end

"""
Evaluate modified Lagrangian.
Should be called before after constraint!(::ElectronicVars)
"""
function eval_L_tilde!( Ham::Hamiltonian, evars::ElectronicVars; kT=1e-3 )

    E_fermi, Entropy = set_occupations!( Ham, kT )
    #println("E_fermi = ", E_fermi)

    Rhoe = calc_rhoe( Ham, evars.ψ )
    update!( Ham, Rhoe )

    Ham.energies = calc_energies( Ham, evars.ψ )
    Ham.energies.mTS = Entropy

    return sum(Ham.energies)
end

"""
Calculate gradient of modified Lagrangian.
Results are stored in `g_evars`.
"""
function grad_eval_L_tilde!(
    Ham::Hamiltonian,
    evars::ElectronicVars,
    g_evars::ElectronicVars;
    kT=1e-3
)

    ψ = evars.ψ
    E_fermi, Entropy = set_occupations!( Ham, kT )

    Rhoe = calc_rhoe( Ham, ψ )
    update!( Ham, Rhoe )

    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt

    g_ψ = g_evars.ψ
    g_η = g_evars.η
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        i = ik + (ispin - 1)*Nkpt
        g_ψ[i], g_η[i] = calc_grad_Haux( Ham, ψ[i], kT, evars.U[i] )
    end

    return

end


"""
Wrapper for `calc_grad_Haux_prec`.
"""
function calc_primary_search_dirs!(
    Ham::Hamiltonian,
    evars::ElectronicVars,
    Δ_evars::ElectronicVars;
    kT=1e-3,
    κ=0.5
)

    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt

    Δ_ψ = Δ_evars.ψ
    Δ_η = Δ_evars.η
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        i = ik + (ispin - 1)*Nkpt
        Δ_ψ[i], Δ_η[i] = calc_grad_Haux_prec( Ham, evars.ψ[i], kT, κ )
    end

    return

end

"""
Calculate preconditioned gradient in η direction
using expression given in PhysRevB.79.241103 (Freysoldt-Boeck-Neugenbauer)
for primary search direction
"""
function calc_grad_Haux_prec(
    Ham::Hamiltonian,
    ψ::Array{ComplexF64,2},
    kT::Float64,
    κ::Float64
)

    ik = Ham.ik
    ispin = Ham.ispin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    ikspin = ik + (ispin - 1)*Nkpt

    # occupation number for this kpoint
    Focc = @view Ham.electrons.Focc[:,ikspin]
    epsilon = @view Ham.electrons.ebands[:,ikspin]

    Ngw     = size(ψ)[1]
    Nstates = size(ψ)[2]

    # gradients
    g_ψ = zeros(ComplexF64, Ngw, Nstates)
    g_η = zeros(ComplexF64, Nstates, Nstates)

    Hψ = op_H( Ham, ψ )

    # subspace Hamiltonian
    Hsub = Hermitian( ψ' * Hψ )

    # gradient for psi (excluding Focc?)
    for ist = 1:Nstates
        g_ψ[:,ist] = Hψ[:,ist]
        for jst = 1:Nstates
            g_ψ[:,ist] = g_ψ[:,ist] - Hsub[jst,ist]*ψ[:,jst]
        end
        g_ψ[:,ist] = Focc[ist]*g_ψ[:,ist]  # FIXME: in the paper there is no `Focc` factor.
    end
    g_ψ = -Kprec( ik, Ham.pw, g_ψ ) # precondition gradient in ψ direction XXX: the sign is negative

    g_η = copy(Hsub)
    # only diagonal elements are different from Hsub
    for ist = 1:Nstates
        g_η[ist,ist] = κ*( Hsub[ist,ist] - epsilon[ist] )
    end

    return g_ψ, g_η  # set the sign ?
end



# using expression given in PhysRevB.79.241103 (Freysoldt-Boeck-Neugenbauer)
function calc_grad_Haux(
    Ham::Hamiltonian,
    ψ::Array{ComplexF64,2},
    kT::Float64,
    U::Matrix{ComplexF64}
)

    ik = Ham.ik
    ispin = Ham.ispin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    ikspin = ik + (ispin - 1)*Nkpt

    # occupation number for this kpoint and spin
    f = copy(Ham.electrons.Focc[:,ikspin])

    if Ham.electrons.Nspin == 1
        f = 0.5*f # for non spin pol
    end
    epsilon = @view Ham.electrons.ebands[:,ikspin]

    Ngw     = size(ψ)[1]
    Nstates = size(ψ)[2]

    # gradients
    g_ψ = zeros(ComplexF64, Ngw, Nstates)
    g_η = zeros(ComplexF64, Nstates, Nstates)

    Hψ = op_H( Ham, ψ )

    # subspace Hamiltonian
    Hsub = Hermitian( ψ' * Hψ )

    # Equation (22)
    # gradient for psi
    for ist = 1:Nstates
        g_ψ[:,ist] = Hψ[:,ist]
        for jst = 1:Nstates
            g_ψ[:,ist] = g_ψ[:,ist] - Hsub[jst,ist]*ψ[:,jst]
        end
        g_ψ[:,ist] = f[ist]*g_ψ[:,ist]
    end

    #println("f       = ", f)
    #println("epsilon = ", epsilon)

    # Equation (24)
    dF_dmu = 0.0
    for ist = 1:Nstates
        dF_dmu = dF_dmu + ( real(Hsub[ist,ist]) - epsilon[ist] ) * f[ist] * (1.0 - f[ist])
    end
    #if Ham.electrons.Nspin == 1
    #    dF_dmu = 2.0*dF_dmu/kT  # XXX multiply by 2 (for non spin pol ?)
    #else
    #    dF_dmu = dF_dmu/kT
    #end
    dF_dmu = dF_dmu/kT
    #@printf("%3d dF_dmu = %18.10f\n", ik, dF_dmu)

    # Equation (19)
    dmu_deta = zeros(Nstates)
    # ss is the denominator of dmu_deta
    ss = 0.0
    for ist = 1:Nstates
        ss = ss + f[ist]*(1.0 - f[ist])
    end
    #@printf("%3d ss dmu_deta = %18.10f\n", ik, ss)

    # including U
    #SMALL = 1e-10
    #if abs(ss) > SMALL
    #    for m = 1:Nstates
    #        ff = 0.0
    #        for n = 1:Nstates
    #            ff = ff + f[n]*(1.0 - f[n]) * real( conj(U[n,m]) * U[m,n] )
    #        end
    #        dmu_deta[m] = ff/ss
    #    end
    #end

    SMALL = 1e-8
    if abs(ss) > SMALL
        for ist = 1:Nstates
            dmu_deta[ist] = f[ist]*(1.0 - f[ist])/ss
        end
    end
    #println("dmu_deta = ", dmu_deta)

    # diagonal of Equation (23)
    for ist = 1:Nstates
        term1 = -( Hsub[ist,ist] - epsilon[ist] ) * f[ist] * ( 1.0 - f[ist] )/kT #* 2 # XXX: non spin pol
        term2 = dmu_deta[ist]*dF_dmu
        #println("")
        #println("ist = ", ist)
        #@printf("Hsub - epsilon = %18.10f\n", real(Hsub[ist,ist]) - epsilon[ist])
        #@printf("term1          = %18.10f\n", real(term1))
        #@printf("term2          = %18.10f\n", real(term2))
        g_η[ist,ist] = term1 + term2
    end

    # off diagonal of Equation (22)
    for ist = 1:Nstates
        for jst = (ist+1):Nstates
            g_η[ist,jst] = Hsub[ist,jst] * (f[ist] - f[jst]) / (epsilon[ist] - epsilon[jst])
            g_η[jst,ist] = g_η[ist,jst]
        end
    end

    return g_ψ, g_η   # set the sign ? (doesn't matter for linmin?)
end
