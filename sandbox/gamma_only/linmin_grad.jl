function linmin_grad!( Ham::Hamiltonian, psiks, g, d, E_orig; αt = 3e-5 )

    psic = zeros_BlochWavefunc(Ham)
    gt = zeros_BlochWavefunc(Ham)

    Nkspin = length(psiks)
    for i in 1:Nkspin
        psic[i] = psiks[i] + αt*d[i]
        ortho_gram_schmidt!( psic[i] )
    end

    calc_grad!( Ham, psic, gt )

    println("psic[1,1] = ", psic[1][1,1])
    println("psic[2,1] = ", psic[1][2,1])

    println("dot(psic,psic) = ", dot_BlochWavefunc(psic,psic))

    println("gt[1,1] = ", gt[1][1,1])
    println("gt[2,1] = ", gt[1][2,1])

    println("dot(g,g) = ", dot_BlochWavefunc(g,g))
    println("dot(gt,gt) = ", dot_BlochWavefunc(gt,gt))

    denum = 2.0*real( dot(g .- gt, d) )
    println("denum = ", denum)    
    if denum != 0.0
        α = abs( αt * 2.0*real( dot(g, d) )/denum )
    else
        α = 0.0
    end

    #m = dot_BlochWavefunc( g, d )
    #do_step!( psic, α, d )
    #E_trial = calc_energies!( Ham, psic )
    #dE = E_orig - E_trial
    #println("α = ", α)

    return true, α

end