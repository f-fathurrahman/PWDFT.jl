function quadratic_line_minimizer!(Ham, psiks, d, E, g, α_t, α)

    Nkspin = length(psiks)
    Δ = zeros(Nkspin)

    # compute directional derivative
    for i in 1:Nkspin
        Δ[i] = 2.0*real( sum( conj(d[i]).*g[i] ) )
        #Δ[i] = 2.0*real( tr(d[i]' * g[i]) )
        #println("Δ[i] = ", Δ[i])
        #if Δ[i] > 0.0
        #    println("Bad step direction!")
        #    α[i] = 0.0
        #end
    end

    println("quadratic_line_minimizer: α_t = ", α_t)
    if sum(Δ) > 0
        println("Wrong direction !!!!")
    end

    psiks_t = psiks + α_t*d
    E_trial = obj_function!( Ham, psiks_t )

    #for i in 1:Nkspin
    #    num = E_trial - ( E + α_t*Δ[i] )
    #    curvature = ( E_trial - ( E + α_t*Δ[i] ) ) /α_t^2
    #    #println("curvature = ", curvature, " num = ", num)
    #    α[i] = -Δ[i]/(2*curvature)
    #    if α[i] < 0
    #        println("Wrong sign of α")
    #        α[i] = α_t # set to zero?
    #    end
    #end

    ss = sum(Δ)
    curvature = ( E_trial - ( E + α_t*ss ) ) /α_t^2
    #println("curvature = ", curvature, " num = ", num)
    cc = -ss/(2*curvature)
    if cc < 0
        println("Wrong sign of α")
        α[:] .= α_t
        α_t = 3.0*α_t
    else
        α[:] .= cc
    end

    return α_t
end
