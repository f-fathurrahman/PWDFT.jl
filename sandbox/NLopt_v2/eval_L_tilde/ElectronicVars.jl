mutable struct ElectronicVars
    ψ::BlochWavefunc
    η::Array{Matrix{ComplexF64},1}
end

function rand_ElectronicVars( Ham::Hamiltonian )

    ψ = rand_BlochWavefunc(Ham)

    Nstates = Ham.electrons.Nstates
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin
    Nkspin = Nkpt*Nspin

    η = Array{Matrix{ComplexF64},1}(undef,Nkspin)
    for i in 1:Nkspin
        η[i] = rand( ComplexF64, Nstates, Nstates )
        η[i] = 0.5*( η[i] + η[i]' )
    end

    return ElectronicVars(ψ, η)
end

"""
Rotate auxiliary Hamiltonian η.
Set Ham.electrons.ebands to eigenvalues of η.
Rotate orthonormalize and rotate ψ.
"""
function constraint!( Ham::Hamiltonian, e::ElectronicVars )

    ψ = e.ψ
    η = e.η

    U = copy(η)
    λ = zeros( Float64, size(η[1],1) ) # eigenvalues
    for i in 1:length(U)
        λ, U[i] = eigen(η[i])
        Ham.electrons.ebands[:,i] = λ
        η[i] = diagm( 0 => λ ) # rotate η
        ortho_sqrt!( ψ[i] )
        ψ[i] = ψ[i] * U[i]
    end

    return
end

function zeros_ElectronicVars( Ham::Hamiltonian )

    ψ = zeros_BlochWavefunc(Ham)

    Nstates = Ham.electrons.Nstates
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin
    Nkspin = Nkpt*Nspin

    η = Array{Matrix{ComplexF64},1}(undef,Nkspin)
    for i in 1:Nkspin
        η[i] = zeros( ComplexF64, Nstates, Nstates )
    end

    return ElectronicVars(ψ, η)
end

import LinearAlgebra: dot
function dot( a::ElectronicVars, b::ElectronicVars )

    N = length(a.ψ)
    v_ψ = zeros(N)
    v_η  = zeros(N)

    for i in 1:N
        v_ψ[i] = real( sum(conj(a.ψ[i]) .* b.ψ[i]) )
        v_η[i]  = real( sum(conj(a.η[i]) .* b.η[i]) )
    end
    return v_ψ, v_η
end


function axpy!( a::Float64, b::Float64, x::ElectronicVars, y::ElectronicVars )
    Nkspin = length(x.ψ)
    # update ψ and η
    for i in 1:Nkspin
        x.ψ[i] = x.ψ[i] + a*y.ψ[i]
        x.η[i] = x.η[i] + b*y.η[i]
        x.η[i] = 0.5*( x.η[i] + x.η[i]' ) # or use previous U_Haux ?
    end
    return
end

function axpy!(
    a::Vector{Float64},
    b::Vector{Float64},
    x::ElectronicVars,
    y::ElectronicVars
)
    Nkspin = length(x.ψ)
    # update ψ and η
    for i in 1:Nkspin
        x.ψ[i] = x.ψ[i] + a[i]*y.ψ[i]
        x.η[i] = x.η[i] + b[i]*y.η[i]
        x.η[i] = 0.5*( x.η[i] + x.η[i]' ) # or use previous U_Haux ?
    end
    return
end

import PWDFT: print_ebands

"""
Wrapper for `print_ebands`.
FIXME: Should be added to `PWDFT.jl`.
"""
function print_ebands( Ham::Hamiltonian )
    print_ebands( Ham.electrons, Ham.pw.gvecw.kpoints, unit="eV" )
end

import Base: copy
function copy( evars::ElectronicVars )
    return ElectronicVars( copy(evars.ψ), copy(evars.η) )
end

"""
Print the first η of an instance of ElectronicVars
"""
function print_Haux( e::ElectronicVars, header::String )
    println()
    println(header)
    println("\nreal part\n")
    display(real(e.η[1]))
    println("\n\nimaginary part\n")
    display(imag(e.η[1]))
    println()
end
