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

function rotate!( e::ElectronicVars )

    psiks = e.psiks
    Haux = e.Haux

    U_Haux = copy(Haux)
    λ = zeros(Float64, size(Haux[1],1) )
    for i in 1:length(U_Haux)
        λ, U_Haux[i] = eigen( Haux[i] )
        Haux[i] = diagm( 0 => λ ) # rotate Haux
        ortho_sqrt!(psiks[i])
        psiks[i] = psiks[i]*U_Haux[i] # rotate psiks
    end

    return
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

import LinearAlgebra: dot
function dot( a::ElectronicVars, b::ElectronicVars )

    N = length(a.psiks)
    v_psiks = zeros(N)
    v_Haux  = zeros(N)

    for i in 1:N
        v_psiks[i] = real( sum(conj(a.psiks[i]) .* b.psiks[i]) )
        v_Haux[i]  = real( sum(conj(a.Haux[i]) .* b.Haux[i]) )
    end
    return v_psiks, v_Haux
end