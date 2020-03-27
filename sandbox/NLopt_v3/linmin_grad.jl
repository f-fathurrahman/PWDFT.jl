function linmin_grad!( Ham, psiks, g, Kg, d; αt = 3e-5 )

    psic = zeros_BlochWavefunc(Ham)
    gt = zeros_BlochWavefunc(Ham)
    Kgt = zeros_BlochWavefunc(Ham)

    Nkspin = length(psiks)
    for i in 1:Nkspin
        psic[i] = psiks[i] + αt*d[i]
        ortho_sqrt!( psic[i] )
    end

    calc_grad!( Ham, psic, gt )

    denum = dot_BlochWavefunc( g .- gt, d )
    if denum != 0.0
        α = abs( αt * dot_BlochWavefunc(g,d)/denum )
    else
        α = 0.0
    end

    println("α = ", α)
    return true, α
#    # Update wavefunction
#    for i in 1:Nkspin
#        psiks[i] = psiks[i] + α * d[i]
#        ortho_sqrt!( psiks[i] )
#    end


#    Etot = calc_energies_grad!( Ham, psiks, g, Kg )

end