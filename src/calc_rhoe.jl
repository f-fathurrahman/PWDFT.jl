include("calc_rhoe_inplace.jl")

function calc_rhoe( Ham::Hamiltonian, psiks::BlochWavefunc )
    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Rhoe = zeros(Float64, Npoints, Nspin)
    calc_rhoe!( Ham, psiks, Rhoe )
    return Rhoe
end


