mutable struct PAWVariables
    lm_fact::Int64    # To converge E_xc integrate up to LM = lm_fact * lm_max
    lm_fact_x::Int64  # As above, for gradient corrected functionals
    xlm::Int64        # Additional factor to add to have a good grad.corr.
    radial_grad_style::Int64 # = 0 or 1, algorithm to use for d/dr
    spheres::Vector{PAWAtomicSphere}
    ddd_paw::Array{Float64,3} # This should be mutable
    total_core_energy::Float64
    E_paw_cmp::Array{Float64,3}
    EHxc_paw::Float64
end


#
# Adapted from paw_init.f90 of QE
#
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

    nhmm = floor(Int64, nhm*(nhm+1)/2)
    ddd_paw = zeros(Float64, nhmm, atoms.Natoms, Nspin)

    # determine constant to constant to be added to total energy
    # to get all-electron energy 
    total_core_energy = 0.0
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species
    for ia in 1:Natoms
        isp = atm2species[ia]
        if pspots[isp].is_paw
            total_core_energy += pspots[isp].paw_data.core_energy
        end
    end
    println("total_core_energy = ", total_core_energy)

    E_paw_cmp = zeros(Float64, Natoms, 2, 2)
    EHxc_paw = 0.0

    return PAWVariables(
        lm_fact, lm_fact_x, xlm, radial_grad_style,
        paw_spheres, ddd_paw, total_core_energy,
        E_paw_cmp, EHxc_paw
    )

end



