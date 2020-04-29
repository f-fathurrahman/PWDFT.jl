function linmin_grad_v1!(
    Ham::Hamiltonian,
    evars_::ElecVars, g::ElecGradient, d::ElecGradient, kT::Float64,
    subrot_; αt = 3e-5
)

    # Probably preallocate these outside linmin_grad?
    evars = deepcopy(evars_)
    subrot = deepcopy(subrot_)
    gt = deepcopy(g)
    Kgt = deepcopy(g)

    do_step!( Ham, Ham, αt, evars, d, subrot )
    Etot = compute!( Ham, evars, gt, Kgt, kT, subrot )

    denum = dot_ElecGradient(g - gt, d)
    if denum != 0.0
        α = abs( αt * dot_ElecGradient(g, d) / denum )
    else
        α = 0.0
    end
    return α

end

function linmin_grad_v2!(
    Ham::Hamiltonian,
    evars_::ElecVars, g::ElecGradient, d::ElecGradient, kT::Float64,
    subrot_; αt = 3e-5
)

    # Probably preallocate these outside linmin_grad?
    evars = deepcopy(evars_)
    subrot = deepcopy(subrot_)
    gt = deepcopy(g)
    Kgt = deepcopy(g)

    do_step!( Ham, αt, evars, d, subrot )
    Etot = compute!( Ham, evars, gt, Kgt, kT, subrot )

    denum, denum_aux = dot_ElecGradient_v2(g - gt, d)
    num, num_aux = dot_ElecGradient_v2(g, d)
    if denum != 0.0
        α = abs( αt * num/denum )
    else
        α = 0.0
    end

    if denum_aux != 0.0
        α_aux = abs( αt * num_aux/denum_aux )
    else
        α_aux = 0.0
    end

    return α, α_aux

    #denum = dot_ElecGradient(g - gt, d)
    #if denum != 0.0
    #    α = abs( αt * dot_ElecGradient(g, d) / denum )
    #else
    #    α = 0.0
    #end
    #return α

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