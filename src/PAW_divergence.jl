function PAW_divergence!(
    ia,
    atoms, pspots, pspotNL,
    F_lm, div_F_lm, lmaxq_in, lmaxq_out
)
    # F_lm has three components: r, θ, ϕ
#=

    !! Compute divergence of a vector field (actually the hamiltonian).  
    !! It is assumed that:  
    !! 1. the input function is multiplied by \(r^2\);  
    !! 2. the output function is multiplied by \(r^2\) too.
    !
=#

    isp = atoms.atm2species[ia]
    Nrmesh = pspots[isp].Nr
    r = pspots[isp].r
    Nspin = size(F_lm, 4)

    ylm = pspotNL.paw.spheres[isp].ylm
    dylmp = pspotNL.paw.spheres[isp].dylmp
    dylmt = pspotNL.paw.spheres[isp].dylmt
    sin_th = pspotNL.paw.spheres[isp].sin_th
    cos_th = pspotNL.paw.spheres[isp].cos_th
    nx = pspotNL.paw.spheres[isp].nx

    div_F_rad = zeros(Float64, Nrmesh, nx, Nspin) # div(F) on radial grid
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
            for lm in 2:lmaxq_in^2
                @views aux[1:Nrmesh] .+= dylmp[ix,lm] .* F_lm[1:Nrmesh,2,lm,ispin]
                # as for PAW_gradient this is already present in dylmp --^
            end
            @views div_F_rad[1:Nrmesh,ix,ispin] .= aux[1:Nrmesh]
        end
    end
    #
    # theta component
    #
    for ispin in 1:Nspin
        for ix in 1:nx
            fill!(aux, 0.0)
            # this derivative has a spherical component too!
            for lm in 1:lmaxq_in^2
                @views aux[1:Nrmesh] .+= F_lm[1:Nrmesh,3,lm,ispin] * (dylmt[ix,lm]*sin_th[ix] + 2.0*ylm[ix,lm]*cos_th[ix])
            end
            @views div_F_rad[1:Nrmesh,ix,ispin] .+= aux[1:Nrmesh] # accumulate
        end
    end
    #
    # Convert what I have done so far to Y_lm
    PAW_rad2lm!( ia, atoms, pspotNL, lmaxq_out, div_F_rad, div_F_lm )
    #
    # Multiply by 1/r**3: 1/r is for theta and phi components only
    # 1/r**2 is common to all the three components.
    for ispin in 1:Nspin
        for lm in 1:lmaxq_out^2
            @views div_F_lm[1:Nrmesh,lm,ispin] .= div_F_lm[1:Nrmesh,lm,ispin] ./ r[1:Nrmesh].^2
        end
    end
    #
    # Compute partial radial derivative d/dr
    for ispin in 1:Nspin
        for lm in 1:lmaxq_out^2
            # Derive along \hat{r} (F already contains a r**2 factor, otherwise
            # it may be better to expand (1/r**2) d(A*r**2)/dr = dA/dr + 2A/r)
            @views radial_gradient_AE!( r, F_lm[1:Nrmesh,1,lm,ispin], aux )
            # Sum it in the divergence: it is already in the right Y_lm form
            @views aux[1:Nrmesh] .= aux[1:Nrmesh] ./ r[1:Nrmesh].^2
            @views div_F_lm[1:Nrmesh,lm,ispin] .= div_F_lm[1:Nrmesh,lm,ispin] .+ aux[1:Nrmesh]
        end
    end

    return
end
