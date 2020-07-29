function linmin_grad!( Ham::HamiltonianGamma,
    psis::BlochWavefuncGamma,
    g::BlochWavefuncGamma,
    d::BlochWavefuncGamma;
    αt = 3e-5
)

    psic = zeros_BlochWavefuncGamma(Ham)
    gt = zeros_BlochWavefuncGamma(Ham)

    Nspin = length(psis)
    for i in 1:Nspin
        psic.data[i] = psis.data[i] + αt*d.data[i]
        #ortho_GS_gamma!(psic.data[i])
        ortho_sqrt_gamma!(psic.data[i])
    end
    
    #println("dot(psic) = ", dot_BlochWavefuncGamma(psic,psic))

    Rhoe = calc_rhoe(Ham, psic)
    update!(Ham, Rhoe)

    for i in 1:Nspin
        Ham.ispin = i
        calc_grad!( Ham, psic.data[i], gt.data[i] )
    end

    denum = 2.0*real( dot_BlochWavefuncGamma(g - gt, d) )
    if denum != 0.0
        α = abs( αt * 2.0*real( dot_BlochWavefuncGamma(g, d) )/denum )
    else
        α = 0.0
    end

    return α

end