function test_shifted_Hamiltonian( σ::Float64, Ham, psiks_old, psiks )

    Hpsi = op_H( Ham, psiks )
    term1 = copy( psiks )  # similar doesn't work

    for i in 1:length(psiks)
        term1[i] = Hpsi[i] - σ*psiks_old[i] * ( psiks_old[i]' * psiks[i] )
    end

    #println("Pass here")

    return
end