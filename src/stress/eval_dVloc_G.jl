#
# dvloc = D Vloc (g^2) / D g^2 = (1/2g) * D Vloc(g) / D g
#
function eval_dVloc_G!(
    psp::PsPot_UPF,
    G2_shells::Vector{Float64},
    dVloc_G::AbstractVector{Float64}
)
    # NOTE: This should be in PsPot_UPF constructor?
    r = psp.r
    Nr_full = psp.Nr
    Nr = Nr_full
    RCUT = 10.0
    for i in 1:Nr_full
        if r[i] > RCUT
            Nr = i
            break
        end 
    end
    # Force Nr to be odd number
    Nr = 2*floor(Int64, (Nr + 1)/2) - 1

    rab = psp.rab
    Vloc_at = psp.V_local
    zval = psp.zval
    Ngl = length(G2_shells)

    fill!(dVloc_G, 0.0)
    aux = zeros(Float64, Nr)
    aux1 = zeros(Float64, Nr)

    # XXX This should always be satisfied for our case
    # The first G-vector is (0,0,0)
    if G2_shells[1] < 1e-8
        # first the G = 0 term
        dVloc_G[1] = 0.0
        igl0 = 2
    else
        igl0 = 1
    end

    # here the G != 0 terms, we first compute the part of the integrand 
    # function independent of |G| in real space
    for ir in 1:Nr
       aux1[ir] = r[ir] * Vloc_at[ir] + zval * erf(r[ir])
    end

    for igl in igl0:Ngl
        Gx = sqrt( G2_shells[igl] )
        # and here we perform the integral, after multiplying for the |G|
        # dependent  part
        #
        # DV(g)/Dg = Integral of r (Dj_0(gr)/Dg) V(r) dr
        for ir in 1:Nr
            aux[ir] = aux1[ir] * ( r[ir] * cos(Gx*r[ir])/Gx - sin(Gx*r[ir])/Gx^2)
        end
        vlcp = PWDFT.integ_simpson(Nr, aux, rab)
        # DV(g^2)/Dg^2 = (DV(g)/Dg)/2g
        vlcp *= 4π / 2.0 / Gx # NOTE: does not include CellVolume factor
        # subtract the long-range term
        G2a = 0.25*G2_shells[igl]
        vlcp += 4π * zval * exp(-G2a) * (G2a + 1.0) / G2_shells[igl]^2
        dVloc_G[igl] = vlcp
    end
    return
end