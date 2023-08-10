function PAW_divergence!(
    i, F_lm, div_F_lm, lmaxq_in, lmaxq_out
)
    # F_lm has three components: r, θ, ϕ


#=

    !! Compute divergence of a vector field (actually the hamiltonian).  
    !! It is assumed that:  
    !! 1. the input function is multiplied by \(r^2\);  
    !! 2. the output function is multiplied by \(r^2\) too.
    !
    USE constants,              ONLY : sqrtpi, fpi, e2
    USE noncollin_module,       ONLY : nspin_gga
    USE lsda_mod,               ONLY : nspin
    USE atom,                   ONLY : g => rgrid
    !
    TYPE(paw_info), INTENT(IN) :: i
    !! atom's minimal info
    INTEGER, INTENT(IN) :: lmaxq_in
    !! max angular momentum to derive (divergence is reliable up to lmaxq_in-2)
    INTEGER, INTENT(IN) :: lmaxq_out
    !! max angular momentum to reconstruct for output
    REAL(DP), INTENT(IN) :: F_lm(i%m,3,lmaxq_in**2,nspin_gga)
    !! Y_lm expansion of F
    REAL(DP), INTENT(OUT):: div_F_lm(i%m,lmaxq_out**2,nspin_gga)
    !! div(F) 
=#

    div_F_rad = zeros(Float64, Nrmesh, nx, Nspin) # div(F) on rad. grid
    aux = zeros(Float64, Nrmesh)
    aux2 = zeros(Float64, Nrmesh)

    #
    # This is the divergence in spherical coordinates:
    #     {1 \over r^2}{\partial ( r^2 A_r ) \over \partial r} 
    #   + {1 \over r\sin\theta}{\partial \over \partial \theta} (  A_\theta\sin\theta )
    #   + {1 \over r\sin\theta}{\partial A_\phi \over \partial \phi}
    #
    # The derivative sum_LM d(Y_LM sin(theta) )/dtheta will be expanded as:
    # sum_LM ( Y_lm cos(theta) + sin(theta) dY_lm/dtheta )
    #
    # The radial component of the divergence is computed last, for practical reasons

    #
    # phi component
    #
    for ispin in 1:Nspin
        for ix in 1:nx
            fill!(aux, 0.0)
            # this derivative has no spherical component, so lm starts from 2
            for lm in 2:(lmaxq_in^2)
                @views aux[1:Nrmesh] .+= dylmp[ix,lm] .* F_lm[1:Nrmesh,2,lm,ispin]
                # as for PAW_gradient this is already present in dylmp --^
            end
            div_F_rad[1:Nrmesh,ix,ispin] .= aux[1:Nrmesh]
        end
    end
    #
    # theta component
    #
    for ispin in 1:Nspin
        for ix in 1:nx
            fill!(aux)
            # this derivative has a spherical component too!
            for lm in 1:(lmaxq_in^2)
                aux[1:Nrmesh] = aux[1:i%m] + F_lm[1:Nrmesh,3,lm,ispin] * (dylmt[ix,lm]*sin_th[ix] +
                                2.0*ylm[ix,lm] * cos_th[ix])
            end
            div_F_rad[1:Nrmesh,ix,is] += aux[1:Nrmesh] # accumulate
        end
    end
    #
    # Convert what I have done so far to Y_lm
    PAW_rad2lm!( i, div_F_rad, div_F_lm, lmaxq_out )
    # Multiply by 1/r**3: 1/r is for theta and phi componente only
    # 1/r**2 is common to all the three components.
    for ispin in 1:Nspin
        for lm in 1:lmaxq_out^2
            div_F_lm[1:Nrmesh,lm,ispin] = div_F_lm[1:i%m,lm,ispin] / r[1:Nrmesh].^2
        end
    end
    #
    # Compute partial radial derivative d/dr
    for ispin in 1:Nspin
        for lm in 1:lmaxq_out^2
            # Derive along \hat{r} (F already contains a r**2 factor, otherwise
            # it may be better to expand (1/r**2) d(A*r**2)/dr = dA/dr + 2A/r)
            radial_gradient_AE!( r, F_lm[1:Nrmesh,1,lm,ispin], aux )
            # Sum it in the divergence: it is already in the right Y_lm form
            aux[1:Nrmesh] = aux[1:Nrmesh] ./ r[1:Nrmesh].^2
            div_F_lm[1:Nrmesh,lm,is] = div_F_lm[1:Nrmesh,lm,is] + aux(1:i%m)
        end
    end

    return
end function