using PWDFT

include("my_trscf.jl")
include("test_shifted_Hamiltonian.jl")

mutable struct ShiftedHamiltonian
    Ham::Ref{Hamiltonian}
    σ::Float64
    psiks_old::BlochWavefunc
end

function op_H( sHam::ShiftedHamiltonian, psi::BlochWavefunc )
end

function test_shifted_Hamiltonian( σ::Float64, Ham, psiks_old, psiks )

    Hpsi = op_H( Ham, psiks )
    term1 = copy( psiks )  # similar doesn't work

    for i in 1:length(psiks)
        term1[i] = Hpsi[i] - σ*psiks_old[i] * ( psiks_old[i]' * psiks[i] )
    end

    return
end

function main()

    Ham = create_Ham_Si_fcc()

    psiks_old = rand_BlochWavefunc( Ham )
    Rhoe = calc_rhoe( Ham, psiks_old )

    update!( Ham, Rhoe )

    psiks = rand_BlochWavefunc( Ham )
    @time test_shifted_Hamiltonian( 0.5, Ham, psiks_old, psiks )
    @time test_shifted_Hamiltonian( 0.5, Ham, psiks_old, psiks )

end

main()