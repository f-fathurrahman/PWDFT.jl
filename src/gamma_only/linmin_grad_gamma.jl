function linmin_grad!(
    Ham::HamiltonianGamma, psis::BlochWavefuncGamma,
    g::BlochWavefuncGamma, d::BlochWavefuncGamma,
    psic::BlochWavefuncGamma, gt::BlochWavefuncGamma;
    αt = 3e-5
)

    Nspin = length(psis)
    for i in 1:Nspin
        psic.data[i] = psis.data[i] + αt*d.data[i]
        ortho_sqrt_gamma!(psic.data[i])
    end

    Rhoe = calc_rhoe(Ham, psic)
    update!(Ham, Rhoe)

    for i in 1:Nspin
        Ham.ispin = i
        calc_grad!( Ham, psic.data[i], gt.data[i] )
    end

    denum = real( dot(g - gt, d) )
    if denum != 0.0
        α = abs( αt * real( dot(g, d) )/denum )
    else
        α = 0.0
    end

    return α

end