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
    Ngw = Ham.pw.gvecw.Ngw
    for ist in 1:Nstates
        for igw in 1:Ngw
            g[igw,ist] = Focc[ist,ispin] * Hψ[igw,ist]
        end
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
    Hsub[:,:] = ψ' * Hψ
    
    v1 = zeros(ComplexF64,Nstates)
    v2 = zeros(ComplexF64,Nstates)
    for ist in 1:Nstates
        v1[ist] = ψ[1,ist]
        v2[ist] = Hψ[1,ist]
    end

    Hsub[:,:] = Hsub + conj(Hsub) - v1*v2'

    Hψ = Hψ - ψ*Hsub
    Ngw = Ham.pw.gvecw.Ngw
    for ist in 1:Nstates
        for igw in 1:Ngw
            g[igw,ist] = Focc[ist,ispin] * Hψ[igw,ist]
        end
    end

    return

end