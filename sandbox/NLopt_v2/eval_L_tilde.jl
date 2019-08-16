using LinearAlgebra
using Printf
using PWDFT
using Random

const DIR_PWDFT = joinpath(dirname(pathof(PWDFT)),"..")
const DIR_PSP = joinpath(DIR_PWDFT, "pseudopotentials", "pade_gth")

include("create_Ham.jl")

function set_occupations!( Ham, kT )

    Ham.electrons.Focc, E_fermi = calc_Focc(
        Ham.electrons.Nelectrons,
        Ham.pw.gvecw.kpoints.wk,
        kT, Ham.electrons.ebands,
        Ham.electrons.Nspin )

    return E_fermi
end

mutable struct ElectronicVars
    psiks::BlochWavefunc
    Haux::Array{Matrix{ComplexF64},1}
end

function rand_ElectronicVars( Ham::Hamiltonian )

    psiks = rand_BlochWavefunc(Ham)
    
    Nstates = Ham.electrons.Nstates
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin
    Nkspin = Nkpt*Nspin

    Haux = Array{Matrix{ComplexF64},1}(undef,Nkspin)
    for i in 1:Nkspin
        Haux[i] = rand( ComplexF64, Nstates, Nstates )
        Haux[i] = 0.5*( Haux[i] + Haux[i]' )
    end

    return ElectronicVars(psiks, Haux)
end

function zeros_ElectronicVars( Ham::Hamiltonian )

    psiks = zeros_BlochWavefunc(Ham)
    
    Nstates = Ham.electrons.Nstates
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin
    Nkspin = Nkpt*Nspin

    Haux = Array{Matrix{ComplexF64},1}(undef,Nkspin)
    for i in 1:Nkspin
        Haux[i] = zeros( ComplexF64, Nstates, Nstates )
    end

    return ElectronicVars(psiks, Haux)
end


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
    Hsub = psi' * Hpsi

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
    Hsub = psi' * Hpsi

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


function test_eval_L_tilde()

    Random.seed!(1234)

    Ham = create_Ham_atom_Pt_smearing()
    evars = rand_ElectronicVars(Ham)

    println("Before eval_L_tilde!")
    display(real(evars.Haux[1]))
    println()
    display(imag(evars.Haux[1]))
    println()

    Etot = eval_L_tilde!(Ham, evars)
    @printf("Etot = %18.10f\n", Etot)

    println("After eval_L_tilde!")
    display(real(evars.Haux[1]))
    println()
    display(imag(evars.Haux[1]))
    println()

end


function test_grad_eval_L_tilde()
    Random.seed!(1234)
    Ham = create_Ham_atom_Pt_smearing()
    evars = rand_ElectronicVars(Ham)
    g_evars = ElectronicVars( copy(evars.psiks), copy(evars.Haux) )

    Etot = eval_L_tilde!(Ham, evars)
    println("Etot = ", Etot)

    grad_eval_L_tilde!(Ham, evars, g_evars)

    Etot = eval_L_tilde!(Ham, evars)
    println("Etot = ", Etot)
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

function test_SD()
    Random.seed!(1234)

    #Ham = create_Ham_atom_Pt_smearing()
    Ham = create_Ham_Al_fcc_smearing()
    evars = rand_ElectronicVars(Ham)

    g_evars = copy(evars)

    Ham.energies.NN = calc_E_NN( Ham.atoms )
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )

    Etot_old = eval_L_tilde!( Ham, evars )

    α_t = 1e-5
    β_t = 1e-1

    for iter = 1:50
        
        grad_eval_L_tilde!( Ham, evars, g_evars )

        axpy!( -α_t, -β_t, evars, g_evars )

        Etot = eval_L_tilde!( Ham, evars )

        @printf("%8d %18.10f %18.10e\n", iter, Etot, Etot_old - Etot)
        #print_ebands( Ham )

        Etot_old = Etot
    end
end
#@time test_SD()

function precond_grad!( Ham, g, Kg; Kscalar=1.0 )
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    for ispin in 1:Nspin, ik in 1:Nkpt
        ikspin = ik + (ispin-1)*Nkpt
        Kg.psiks[ikspin] = Kprec( ik, Ham.pw, g.psiks[ikspin] )
        Kg.Haux[ikspin] = Kscalar*g.Haux[ikspin]
    end
    return
end

function calc_beta_CG!( g, g_old, Kg, Kg_old, β, β_Haux )

    Nkspin = length(g.psiks)

    for i in 1:Nkspin
        ss = real(sum(conj(g_old.psiks[i]).*Kg_old.psiks[i]))
        if abs(ss) >= 1e-10
            β[i] = real(sum(conj(g.psiks[i]-g_old.psiks[i]).*Kg.psiks[i]))/ss
        else
            β[i] = 0.0
        end
        if β[i] < 0.0
            β[i] = 0.0
        end
        ss = real(sum(conj(g_old.Haux[i]).*Kg_old.Haux[i]))
        if abs(ss) >= 1e-10
            β_Haux[i] = real(sum(conj(g.Haux[i]-g_old.Haux[i]).*Kg.Haux[i]))/ss
        else
            β_Haux[i] = 0.0
        end
        if β_Haux[i] < 0.0
            β_Haux[i] = 0.0
        end
    end
    return
end


function calc_search_dirs!( d, Kg, d_old, β, β_Haux )
    Nkspin = length(d.psiks)
    for i in 1:Nkspin
        d.psiks[i] = Kg.psiks[i] + β[i] * d_old.psiks[i]
        d.Haux[i] = Kg.Haux[i] + β_Haux[i] * d_old.Haux[i]
    end
    return
end

function trial_evars!( ec, e, d, α_t, α_t_aux )
    Nkspin = length(e.psiks)
    for i in 1:Nkspin
        ec.psiks[i] = e.psiks[i] + α_t*d.psiks[i]
        ec.Haux[i] = e.Haux[i] + α_t_aux*d.Haux[i]
        ec.Haux[i] = 0.5*( ec.Haux[i] + ec.Haux[i]' )
    end
    return
end

function calc_alpha_CG!(
    g::ElectronicVars, gt::ElectronicVars, d::ElectronicVars,
    α_t::Float64, α_t_aux::Float64,
    α, α_aux
)

    Nkspin = length(g.psiks)

    for i in 1:Nkspin
        
        denum = real(sum(conj(g.psiks[i]-gt.psiks[i]).*d.psiks[i]))
        #if abs(denum) <= 1e-6
        if denum != 0.0
            α[i] = abs( α_t*real(sum(conj(g.psiks[i]).*d.psiks[i]))/denum )
        else
            α[i] = 0.0
        end

        denum_aux = real(sum(conj(g.Haux[i]-gt.Haux[i]).*d.Haux[i]))
        #if abs(denum) <= 1e-6
        if denum_aux != 0.0
            α_aux[i] = abs( α_t_aux*real(sum(conj(g.Haux[i]).*d.Haux[i]))/denum_aux )
        else
            α_aux[i] = 0.0
        end

    end
    return
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

function test_CG()
    Random.seed!(1234)

    #Ham = create_Ham_atom_Pt_smearing()
    Ham = create_Ham_Pt_fcc_smearing()
    #Ham = create_Ham_Al_fcc_smearing()
    evars = rand_ElectronicVars(Ham)

    evarsc = copy(evars)
    g_evars = copy(evars)
    Kg_evars = copy(evars)
    gt_evars = copy(evars)
    g_old_evars = copy(evars)
    Kg_old_evars = copy(evars)
    d_evars = copy(evars)
    d_old_evars = zeros_ElectronicVars(Ham)

    Nkspin = length(evars.psiks)

    β = zeros(Float64,Nkspin)
    β_Haux = zeros(Float64,Nkspin)

    α = zeros(Float64,Nkspin)
    α_Haux = zeros(Float64,Nkspin)

    Ham.energies.NN = calc_E_NN( Ham.atoms )
    Ham.energies.PspCore = calc_PspCore_ene( Ham.atoms, Ham.pspots )

    Etot_old = eval_L_tilde!( Ham, evars )

    println("Etot_old = ", Etot_old)

    α_t = 1e-5

    for iter = 1:50
        
        grad_eval_L_tilde!( Ham, evars, g_evars )
        #print_Haux( evars, "evars after grad_eval_L_tilde")
        #print_Haux( g_evars, "g_evars after grad_eval_L_tilde")

        calc_primary_search_dirs!( Ham, evars, Kg_evars )
        #print_Haux( Kg_evars, "Kg_evars after grad_eval_L_tilde")

        if iter > 1
            calc_beta_CG!( g_evars, g_old_evars, Kg_evars, Kg_old_evars, β, β_Haux )
        end
        println("β = ", β)
        println("β_Haux = ", β_Haux)

        calc_search_dirs!( d_evars, Kg_evars, d_old_evars, β, β_Haux )
        #print_Haux( d_evars, "d_evars after calc_search_dirs!")

        trial_evars!( evarsc, evars, d_evars, α_t, α_t )
        #print_Haux( evarsc, "evarsc after trial_evars!")

        #print_Haux( evarsc, "evarsc before grad_eval_L_tilde")
        grad_eval_L_tilde!( Ham, evarsc, gt_evars )
        #print_Haux( evarsc, "evarsc after grad_eval_L_tilde")

        calc_alpha_CG!( evars, gt_evars, d_evars, α_t, α_t, α, α_Haux )

        println("α      = ", α)
        println("α_Haux = ", α_Haux)

        # update evars
        axpy!( α, α_Haux, evars, d_evars )

        #print_Haux( evars, "evars before eval_L_tilde")
        Etot = eval_L_tilde!( Ham, evars )
        #print_Haux( evars, "evars after eval_L_tilde")

        @printf("%8d %18.10f %18.10e\n", iter, Etot, Etot_old - Etot)
        #print_ebands( Ham )

        Etot_old = Etot
        d_old_evars = copy(d_evars)
        g_old_evars = copy(g_evars)
        Kg_old_evars = copy(Kg_evars)

    end

    print_ebands( Ham )


end
@time test_CG()