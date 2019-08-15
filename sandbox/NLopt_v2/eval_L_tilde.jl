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

function axpy!(a, x::ElectronicVars, y::ElectronicVars )
    
    # update psiks
    x.psiks = x.psiks + a*y.psiks
    
    # update Haux
    for i in length(x.Haux)
        x.Haux[i] = x.Haux[i] + a*y.Haux[i]
        x.Haux[i] = 0.5*( x.Haux[i] + x.Haux[i]' ) # or use previous U_Haux
    end

    return
end

function test_SD()
    Random.seed!(1234)

    Ham = create_Ham_atom_Pt_smearing()
    evars = rand_ElectronicVars(Ham)

    g_evars = ElectronicVars( copy(evars.psiks), copy(evars.Haux) )

    Etot_old = eval_L_tilde!(Ham, evars)

    α_t = 1e-5

    for iter = 1:5
        
        grad_eval_L_tilde!( Ham, evars, g_evars )

        axpy!( -α_t, evars, g_evars )

        Etot = eval_L_tilde!( Ham, evars )

        @printf("%8d %18.10f %18.10e\n", iter, Etot, Etot_old - Etot)

        Etot_old = Etot
    end
end
test_SD()
