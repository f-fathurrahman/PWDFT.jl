function calc_energies_grad!(
    Ham::HamiltonianGamma,
    psis::BlochWavefuncGamma,
    g::BlochWavefuncGamma,
    Kg::BlochWavefuncGamma,
    Hsub::Vector{Matrix{ComplexF64}}
)

    Rhoe = calc_rhoe( Ham, psis )
    update!( Ham, Rhoe )

    Ham.energies = calc_energies( Ham, psis )

    Nspin = Ham.electrons.Nspin
    for ispin in 1:Nspin
        Ham.ispin = ispin
        calc_grad!( Ham, psis.data[ispin], g.data[ispin], Hsub[ispin] )
        Kprec!( Ham.pw, g.data[ispin], Kg.data[ispin] )
    end

    return sum( Ham.energies )
end


# Modify Ham and psiks
function calc_energies_only!( Ham, psiks )
    Rhoe = calc_rhoe( Ham, psiks )
    update!( Ham, Rhoe )
    Ham.energies = calc_energies( Ham, psiks )
    return sum( Ham.energies )
end

function Kprec!( pw::PWGridGamma, psi::Array{ComplexF64,2}, Kpsi::Array{ComplexF64,2} )
    Ngw  = size(psi,1)
    Nstates = size(psi,2)
    idx_gw2g = pw.gvecw.idx_gw2g
    G2 = pw.gvec.G2
    for ist = 1:Nstates
        for igw = 1:Ngw
            ig = idx_gw2g[igw]
            Kpsi[igw,ist] = psi[igw,ist] / ( 1.0 + G2[ig] )
        end
    end
    return
end

function constrain_search_dir!( d::BlochWavefuncGamma, psis::BlochWavefuncGamma )
    Nspin = length(psis)
    for i in 1:Nspin
        constrain_search_dir_gamma!( d.data[i], psis.data[i] )
    end
    return
end

function constrain_search_dir_gamma!( d::Array{ComplexF64,2}, psi::Array{ComplexF64,2} )
    C = psi' * d
    Nstates = size(d,2)
    v1g = zeros(ComplexF64,Nstates)
    v2g = zeros(ComplexF64,Nstates)
    for ist in 1:Nstates
        v1g[ist] = psi[1,ist] # psi is the array that is conj transposed in the orig expression
        v2g[ist] = d[1,ist]  # v2g is the array that will be conj transposed
    end
    d[:] = d - psi * ( C + conj(C) - v1g*v2g' )
    return
end

function do_step!( α::Float64, psis::BlochWavefuncGamma, d::BlochWavefuncGamma )
    Nspin = length(psis)
    for i in 1:Nspin
        psis.data[i] = psis.data[i] + α*d.data[i]
        #ortho_GS_gamma!( psis.data[i] )
        ortho_sqrt_gamma!( psis.data[i] )
    end
    return
end
