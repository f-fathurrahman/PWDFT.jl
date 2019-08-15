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

import Base: +
function +( x::ElectronicVars, y::ElectronicVars )
    z = ElectronicVars( copy(x.psiks), copy(y.Haux) )
    for i in 1:length(x.psiks)
        z.psiks[i] = x.psiks[i] + y.psiks[i]
        z.Haux[i]  = x.Haux[i] + y.Haux[i]
        z.Haux[i]  = 0.5*( z.Haux[i] + z.Haux[i]' )
    end
    return z
end

function eval_L_tilde!( Ham::Hamiltonian, evars::ElectronicVars; kT=1e-3 )

    psiks = evars.psiks
    Haux = evars.Haux

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

    print_ebands( Ham.electrons, Ham.pw.gvecw.kpoints )

    return sum(energies)
end

function test_main()

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

    println("After  eval_L_tilde!")
    display(real(evars.Haux[1]))
    println()
    display(imag(evars.Haux[1]))
    println()

end

test_main()