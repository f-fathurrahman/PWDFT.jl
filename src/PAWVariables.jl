struct PAWVariables
    lm_fact::Int64    # To converge E_xc integrate up to LM = lm_fact * lm_max
    lm_fact_x::Int64  # As above, for gradient corrected functionals
    xlm::Int64        # Additional factor to add to have a good grad.corr.
    radial_grad_style::Int64 # = 0 or 1, algorithm to use for d/dr
    spheres::Vector{PAWAtomicSphere}
    ddd_paw::Array{Float64,3}
end

function PAWVariables(atoms::Atoms, pspots, nhm::Int64; is_gga=false, Nspin=1)
    lm_fact = 3
    lm_fact_x = 3
    xlm = 2
    radial_grad_style = 0

    lmax_safe = 0
    lmax_add = 0

    Nspecies = atoms.Nspecies
    paw_spheres = Vector{PAWAtomicSphere}(undef,Nspecies)

    for isp in 1:Nspecies
        psp = pspots[isp]
        if !psp.is_paw
            continue
        end
        if psp.lmax_rho == 0
            # no need for more than one direction, when it is spherical!
            lmax_safe = 0
            lmax_add  = 0
        else
            if is_gga
                # Integrate up to a higher maximum lm if using gradient
                # correction check expression for d(y_lm)/d\theta for details
                lmax_safe = lm_fact_x * psp.lmax_rho
                lmax_add  = xlm
            else
                # no gradient correction:
                lmax_safe = lm_fact * psp.lmax_rho
                lmax_add  = 0 
            end
        end

        println("lmax_safe = ", lmax_safe)
        println("lmax_add  = ", lmax_add)

        paw_spheres[isp] = PAWAtomicSphere(lmax_safe, lmax_add)

    end

    nhmm = floor(Int64, nhm*(nhm+1))
    ddd_paw = zeros(Float64, nhmm, atoms.Natoms, Nspin)

    return PAWVariables(
        lm_fact, lm_fact_x, xlm, radial_grad_style,
        paw_spheres, ddd_paw
    )

end



