import PWDFT: calc_grad

function calc_grad( Ham::HamiltonianGamma, ψ::Array{ComplexF64,2} )
    g = zeros(ComplexF64,size(ψ))
    Nstates = size(ψ,2)
    Hsub = zeros(ComplexF64,Nstates,Nstates)
    calc_grad!( Ham, ψ, g, Hsub )
    return g
end

# No Hsub
function calc_grad!( Ham::HamiltonianGamma, ψ::Array{ComplexF64,2}, g::Array{ComplexF64,2} )

    ispin = Ham.ispin

    Nstates = size(ψ,2)
    Focc = Ham.electrons.Focc

    Hψ = op_H( Ham, ψ )
    Hsub = ψ' * Hψ

    v1 = zeros(ComplexF64,Nstates)
    v2 = zeros(ComplexF64,Nstates)
    for ist in 1:Nstates
        v1[ist] = ψ[1,ist]
        v2[ist] = Hψ[1,ist]
    end

    Hsub = Hsub + conj(Hsub) - v1*v2'
    Hψ = Hψ - ψ*Hsub

    for ist in 1:Nstates
        g[:,ist] = Focc[ist,ispin] * Hψ[:,ist]
    end

    return

end


function calc_grad!(
    Ham::HamiltonianGamma,
    ψ::Array{ComplexF64,2},
    g::Array{ComplexF64,2},
    Hsub::Matrix{ComplexF64}
)

    ispin = Ham.ispin

    Nstates = size(ψ,2)
    Focc = Ham.electrons.Focc

    Hψ = op_H( Ham, ψ )

    #println("Hψ[1,1] = ", Hψ[1,1])

    Hsub[:] = ψ' * Hψ
    
    v1 = zeros(ComplexF64,Nstates)
    v2 = zeros(ComplexF64,Nstates)
    for ist in 1:Nstates
        v1[ist] = ψ[1,ist]
        v2[ist] = Hψ[1,ist]
    end

    Hsub[:] = Hsub + conj(Hsub) - v1*v2'

    #println("Cross term v1*v2'")
    #display(v1*v2'); println()

    Hψ = Hψ - ψ*Hsub

    for ist in 1:Nstates
        g[:,ist] = Focc[ist,ispin] * Hψ[:,ist]
    end

    return

end