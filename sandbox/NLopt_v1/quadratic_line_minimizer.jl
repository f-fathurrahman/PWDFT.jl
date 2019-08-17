function quadratic_line_minimizer!(Ham, psiks, d, E, g, α_t, α)
    
    Nkspin = length(psiks)
    Δ = zeros(Nkspin)
    
    # compute directional derivative
    for i in 1:Nkspin
        Δ[i] = 2.0*real( sum( conj(d[i]).*g[i] ) )
        println("Δ[i] = ", Δ[i])
    end
    
    psiks_t = psiks + α_t*d
    E_trial = obj_function!( Ham, psiks_t )
    
    for i in 1:Nkspin
        num = E_trial - ( E + α_t*Δ[i] )
        curvature = ( E_trial - ( E + α_t*Δ[i] ) ) /α_t^2
        println("curvature = ", curvature, " num = ", num)
        α[i] = -Δ[i]/(2*curvature)
    end
    return
end