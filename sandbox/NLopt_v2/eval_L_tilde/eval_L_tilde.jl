using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("../create_Ham.jl")
include("subspace_rotation.jl")
include("ElectronicVars.jl")

function eval_L_tilde!( Ham::Hamiltonian, evars::ElectronicVars; kT=1e-3, skip_ortho=false )

    psiks = evars.psiks
    Haux = evars.Haux

    if !skip_ortho
        for i in length(psiks)
            ortho_sqrt!(psiks[i])
        end
    end

    U_Haux = copy(Haux)
    for i in 1:length(U_Haux)
        Ham.electrons.ebands[:,i], U_Haux[i] = eigen( Haux[i] )
        Haux[i] = diagm( 0 => Ham.electrons.ebands[:,i] ) # rotate Haux
        psiks[i] = psiks[i]*U_Haux[i] # rotate psiks
    end

    E_fermi = set_occupations!( Ham, kT )
    Entropy = calc_entropy(
        Ham.pw.gvecw.kpoints.wk,
        kT,
        Ham.electrons.ebands,
        E_fermi,
        Ham.electrons.Nspin
    )

    Rhoe = calc_rhoe( Ham, psiks )
    update!( Ham, Rhoe )

    energies = calc_energies( Ham, psiks )
    energies.mTS = Entropy

    return sum(energies)
end


function eval_L_tilde_and_rotate_dirs!(
    Ham::Hamiltonian,
    evars::ElectronicVars,
    d::ElectronicVars;
    kT=1e-3,
    skip_ortho=false
)

    psiks = evars.psiks
    Haux = evars.Haux

    if !skip_ortho
        for i in length(psiks)
            ortho_sqrt!(psiks[i])
        end
    end

    U_Haux = copy(Haux)
    d_U_Haux = copy(Haux)

    for i in 1:length(U_Haux)
        Ham.electrons.ebands[:,i], U_Haux[i] = eigen( Haux[i] )
        Haux[i] = diagm( 0 => Ham.electrons.ebands[:,i] ) # rotate Haux
        psiks[i] = psiks[i]*U_Haux[i] # rotate psiks

        # also rotate search direction
        λ, d_U_Haux[i] = eigen( d.Haux[i] )
        d.Haux[i] = diagm( 0 => λ )
        #ortho_sqrt!(d.psiks[i])  # need this?
        d.psiks[i] = d.psiks[i]*d_U_Haux[i]
        
        # using current U_Haux
        #d.psiks[i] = d.psiks[i]*U_Haux[i]
        #d.Haux[i] = U_Haux[i]' * d.Haux[i] * U_Haux[i]
    end

    E_fermi = set_occupations!( Ham, kT )
    Entropy = calc_entropy(
        Ham.pw.gvecw.kpoints.wk,
        kT,
        Ham.electrons.ebands,
        E_fermi,
        Ham.electrons.Nspin
    )

    Rhoe = calc_rhoe( Ham, psiks )
    update!( Ham, Rhoe )

    energies = calc_energies( Ham, psiks )
    energies.mTS = Entropy

    return sum(energies)
end



function grad_eval_L_tilde!(
    Ham::Hamiltonian,
    evars::ElectronicVars,
    g_evars::ElectronicVars;
    kT=1e-3,
    skip_ortho=false
)

    psiks = evars.psiks
    Haux = evars.Haux

    if !skip_ortho
        for i in length(psiks)
            ortho_sqrt!(psiks[i])
        end
    end

    U_Haux = copy(Haux)
    for i in 1:length(U_Haux)
        Ham.electrons.ebands[:,i], U_Haux[i] = eigen( Haux[i] )
        Haux[i] = diagm( 0 => Ham.electrons.ebands[:,i] ) # rotate Haux
        psiks[i] = psiks[i]*U_Haux[i] # rotate psiks
    end

    E_fermi = set_occupations!( Ham, kT )
    Entropy = calc_entropy(
        Ham.pw.gvecw.kpoints.wk,
        kT,
        Ham.electrons.ebands,
        E_fermi,
        Ham.electrons.Nspin
    )

    Rhoe = calc_rhoe( Ham, psiks )
    update!( Ham, Rhoe )

    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt

    g = g_evars.psiks
    g_Haux = g_evars.Haux
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        i = ik + (ispin - 1)*Nkpt
        g[i], g_Haux[i] = calc_grad_Haux(Ham, psiks[i], kT)
    end

    return

end


function calc_primary_search_dirs!(
    Ham::Hamiltonian,
    evars::ElectronicVars,
    Δ_evars::ElectronicVars;
    kT=1e-3,
    skip_ortho=false,
    κ=0.5
)

    psiks = evars.psiks
    Haux = evars.Haux

    if !skip_ortho
        for i in length(psiks)
            ortho_sqrt!(psiks[i])
        end
    end

    U_Haux = copy(Haux)
    for i in 1:length(U_Haux)
        Ham.electrons.ebands[:,i], U_Haux[i] = eigen( Haux[i] )
        Haux[i] = diagm( 0 => Ham.electrons.ebands[:,i] ) # rotate Haux
        psiks[i] = psiks[i]*U_Haux[i] # rotate psiks
    end

    E_fermi = set_occupations!( Ham, kT )
    Entropy = calc_entropy(
        Ham.pw.gvecw.kpoints.wk,
        kT,
        Ham.electrons.ebands,
        E_fermi,
        Ham.electrons.Nspin
    )

    Rhoe = calc_rhoe( Ham, psiks )
    update!( Ham, Rhoe )

    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt

    Δ = Δ_evars.psiks
    Δ_Haux = Δ_evars.Haux
    for ispin in 1:Nspin, ik in 1:Nkpt
        Ham.ispin = ispin
        Ham.ik = ik
        i = ik + (ispin - 1)*Nkpt
        Δ[i], Δ_Haux[i] = calc_grad_Haux_prec(Ham, psiks[i], kT, κ)
        Δ[i] = -Kprec( ik, Ham.pw, Δ[i] )
    end

    return

end


# using expression given in PhysRevB.79.241103 (Freysoldt-Boeck-Neugenbauer)
# for primary search direction
function calc_grad_Haux_prec(
    Ham::Hamiltonian,
    psi::Array{ComplexF64,2},
    kT::Float64,
    κ::Float64
)

    ik = Ham.ik
    ispin = Ham.ispin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    ikspin = ik + (ispin - 1)*Nkpt

    # occupation number for this kpoint
    epsilon = @view Ham.electrons.ebands[:,ikspin]
    
    Ngw     = size(psi)[1]
    Nstates = size(psi)[2]

    # gradients
    g_psi = zeros(ComplexF64, Ngw, Nstates)
    g_Haux = zeros(ComplexF64, Nstates, Nstates)

    Hpsi = op_H( Ham, psi )

    # subspace Hamiltonian
    Hsub = Hermitian( psi' * Hpsi )

    # gradient for psi (excluding Focc)
    for ist = 1:Nstates
        g_psi[:,ist] = Hpsi[:,ist]
        for jst = 1:Nstates
            g_psi[:,ist] = g_psi[:,ist] - Hsub[jst,ist]*psi[:,jst]
        end
    end

    g_Haux = copy(Hsub)
    # diagonal
    for ist = 1:Nstates
        g_Haux[ist,ist] = κ*( Hsub[ist,ist] - epsilon[ist] )
    end

    return g_psi, g_Haux
end



# using expression given in PhysRevB.79.241103 (Freysoldt-Boeck-Neugenbauer)
function calc_grad_Haux(
    Ham::Hamiltonian,
    psi::Array{ComplexF64,2},
    kT::Float64
)

    ik = Ham.ik
    ispin = Ham.ispin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    ikspin = ik + (ispin - 1)*Nkpt

    # occupation number for this kpoint
    f = @view Ham.electrons.Focc[:,ikspin]
    epsilon = @view Ham.electrons.ebands[:,ikspin]
    
    Ngw     = size(psi)[1]
    Nstates = size(psi)[2]

    # gradients
    g_psi = zeros(ComplexF64, Ngw, Nstates)
    g_Haux = zeros(ComplexF64, Nstates, Nstates)

    Hpsi = op_H( Ham, psi )

    # subspace Hamiltonian
    Hsub = Hermitian( psi' * Hpsi )

    # gradient for psi
    for ist = 1:Nstates
        g_psi[:,ist] = Hpsi[:,ist]
        for jst = 1:Nstates
            g_psi[:,ist] = g_psi[:,ist] - Hsub[jst,ist]*psi[:,jst]
        end
        g_psi[:,ist] = f[ist]*g_psi[:,ist]
    end


    dF_dmu = 0.0
    for ist = 1:Nstates
        dF_dmu = dF_dmu + ( Hsub[ist,ist] - epsilon[ist] ) * 0.5*f[ist] * (1.0 - 0.5*f[ist])
    end
    dF_dmu = dF_dmu/kT

    dmu_deta = zeros(Nstates)
    # ss is the denominator of dmu_deta
    ss = 0.0
    for ist = 1:Nstates
        ss = ss + 0.5*f[ist]*(1.0 - 0.5*f[ist])
    end
    SMALL = 1e-8
    if abs(ss) > SMALL
        for ist = 1:Nstates
            dmu_deta[ist] = 0.5*f[ist]*(1.0 - 0.5*f[ist])/ss
        end
    end

    # diagonal
    for ist = 1:Nstates
        term1 = -( Hsub[ist,ist] - epsilon[ist] ) * 0.5*f[ist] * ( 1.0 - 0.5*f[ist] )/kT
        term2 = dmu_deta[ist]*dF_dmu
        g_Haux[ist,ist] = term1 + term2
    end

    # off diagonal
    for ist = 1:Nstates
        for jst = (ist+1):Nstates
            g_Haux[ist,jst] = Hsub[ist,jst] * (f[ist] - f[jst]) / (epsilon[ist] - epsilon[jst])
            g_Haux[jst,ist] = g_Haux[ist,jst]
        end
    end

    return g_psi, g_Haux
end

function axpy!(a::Float64, b::Float64, x::ElectronicVars, y::ElectronicVars )
    
    Nkspin = length(x.psiks)
    # update psiks and Haux
    for i in 1:Nkspin
        x.psiks[i] = x.psiks[i] + a*y.psiks[i]
        x.Haux[i] = x.Haux[i] + b*y.Haux[i]
        x.Haux[i] = 0.5*( x.Haux[i] + x.Haux[i]' ) # or use previous U_Haux ?
    end

    return
end

function axpy!(
    a::Vector{Float64},
    b::Vector{Float64},
    x::ElectronicVars,
    y::ElectronicVars
)
    
    Nkspin = length(x.psiks)
    # update psiks and Haux
    for i in 1:Nkspin
        x.psiks[i] = x.psiks[i] + a[i]*y.psiks[i]
        x.Haux[i] = x.Haux[i] + b[i]*y.Haux[i]
        x.Haux[i] = 0.5*( x.Haux[i] + x.Haux[i]' ) # or use previous U_Haux ?
    end

    return
end

import PWDFT: print_ebands
function print_ebands( Ham::Hamiltonian )
    print_ebands( Ham.electrons, Ham.pw.gvecw.kpoints )
end

import Base: copy
function copy( evars::ElectronicVars )
    return ElectronicVars( copy(evars.psiks), copy(evars.Haux) )
end

# print the first Haux of an instance of ElectronicVars
function print_Haux( e::ElectronicVars, header::String )
    println()
    println(header)
    println("\nreal part\n")
    display(real(e.Haux[1]))
    println("\n\nimaginary part\n")
    display(imag(e.Haux[1]))
    println()
end

