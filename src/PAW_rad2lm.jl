function PAW_rad2lm!(
    ia,
    atoms, pspotNL,
    lmax_loc::Int64,
    F_rad, F_lm
)
    # Computes:
    # \[ F_{lm}(r) = \int d \Omega\ F(r,\text{th},\text{ph})\ Y_{lm}(\text{th},
    # \text{ph}) \]

    Nspin = size(F_rad, 3)
    Nrmesh = size(F_rad, 1)
    isp = atoms.atm2species[ia]
    sphere = pspotNL.paw.spheres[isp]

    for ispin in 1:Nspin, lm in 1:lmax_loc^2
        F_lm[:,lm,ispin] .= 0.0
        for ix in 1:sphere.nx, j in 1:Nrmesh
             F_lm[j,lm,ispin] += F_rad[j,ix,ispin] * sphere.wwylm[ix,lm]
        end
    end

    return
end