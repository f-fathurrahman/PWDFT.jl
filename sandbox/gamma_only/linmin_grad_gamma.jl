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
        ortho_GS_gamma!(psic.data[i])
    end

    Rhoe = calc_rhoe(Ham, psic)
    update!(Ham, Rhoe)

    for i in 1:Nspin
        Ham.ispin = i
        calc_grad!( Ham, psic.data[i], gt.data[i] )
    end

    #ortho_check( psic )

    #println("psic[1,1] = ", psic.data[1][1,1])
    #println("psic[2,1] = ", psic.data[1][2,1])
    
    #println("dot(psic,psic) = ", 2*real(dot_BlochWavefuncGamma(psic,psic)))

    #println("gt[1,1] = ", gt.data[1][1,1])
    #println("gt[2,1] = ", gt.data[1][2,1])

    #println("dot(g,g) = ", 2*real(dot_BlochWavefuncGamma(g,g)))
    #println("dot(gt,gt) = ", 2*real(dot_BlochWavefuncGamma(gt,gt)))

    denum = 2.0*real( dot_BlochWavefuncGamma(g - gt, d) )
    #println("denum = ", denum)
    if denum != 0.0
        α = abs( αt * 2.0*real( dot_BlochWavefuncGamma(g, d) )/denum )
    else
        α = 0.0
    end

    return α

end