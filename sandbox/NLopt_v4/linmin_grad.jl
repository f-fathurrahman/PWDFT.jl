function linmin_grad!(
    Ham::Hamiltonian,
    evars_::ElecVars, g::ElecGradient, d::ElecGradient, kT::Float64,
    rotPrev_, rotPrevC_, rotPrevCinv_; αt = 3e-5
)

    # Probably preallocate these outside linmin_grad?
    evars = deepcopy(evars_)
    rotPrev = deepcopy(rotPrev_)
    rotPrevC = deepcopy(rotPrevC_)
    rotPrevCinv = deepcopy(rotPrevCinv_)
    gt = deepcopy(g)
    Kgt = deepcopy(g)

    do_step!( αt, evars, d, rotPrev, rotPrevC, rotPrevCinv )
    Etot = compute!( Ham, evars, gt, Kgt, kT, rotPrevCinv, rotPrev )

    denum = dot_ElecGradient(g - gt, d)
    if denum != 0.0
        α = abs( αt * dot_ElecGradient(g, d) / denum )
    else
        α = 0.0
    end

    return α

end

function linmin_grad!( Ham, psiks::BlochWavefunc, g::BlochWavefunc, d::BlochWavefunc; αt = 3e-5 )

    psic = zeros_BlochWavefunc(Ham)
    gt = zeros_BlochWavefunc(Ham)

    Nkspin = length(psiks)
    for i in 1:Nkspin
        psic[i] = psiks[i] + αt*d[i]
        ortho_sqrt!( psic[i] )
    end

    calc_grad!( Ham, psic, gt )

    #denum = dot_BlochWavefunc( g .- gt, d )
    denum = 2.0*real( dot(g .- gt, d) )
    if denum != 0.0
        α = abs( αt * 2.0*real( dot(g, d) )/denum )
    else
        α = 0.0
    end

    return α

end