function linmin_grad!( Ham::Hamiltonian, psiks, g, d, E_orig; αt = 3e-5 )

    psic = zeros_BlochWavefunc(Ham)
    gt = zeros_BlochWavefunc(Ham)

    Nkspin = length(psiks)
    for i in 1:Nkspin
        psic[i] = psiks[i] + αt*d[i]
        ortho_gram_schmidt!( psic[i] )
    end

    calc_grad!( Ham, psic, gt )

    denum = 2.0*real( dot(g .- gt, d) )
    if denum != 0.0
        α = abs( αt * 2.0*real( dot(g, d) )/denum )
    else
        α = 0.0
    end

    return α

end