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
    Hsub = Hsub + conj(Hsub)
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
    Hsub[:] = ψ' * Hψ
    Hsub = Hsub + conj(Hsub)
    Hψ = Hψ - ψ*Hsub

    for ist in 1:Nstates
        g[:,ist] = Focc[ist,ispin] * Hψ[:,ist]
    end

    return

end