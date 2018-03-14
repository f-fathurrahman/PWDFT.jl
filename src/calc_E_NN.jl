#
# Copyright (C) 2001 PWSCF group
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#
# modified by Fadjar Fathurrahman
#

using SpecialFunctions: erfc

include("rgen.jl")

"""
Calculates Ewald energy with both G- and R-space terms.
Determines optimal alpha. Should hopefully work for any structure.
"""
function calc_E_NN( pw::PWGrid, atoms::Atoms, Zv::Array{Float64,1} )

    Sf = calc_strfact( atoms, pw )

    atpos = atoms.positions
    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species

    Ns = pw.Ns
    LL = pw.LatVecs
    Ng = pw.gvec.Ng
    G2 = pw.gvec.G2
    omega = pw.Ω

    ewaldg = 0.0
    ewaldsr = 0.0

    dtau = zeros(Float64,3)

    # REAL(DP) :: charge, ewaldg, ewaldr, dtau (3), alpha, &
    # r (3, mxr), r2 (mxr), rmax, rr, upperbound, fact

    at = pw.LatVecs'  # assuming transpose
    bg = pw.RecVecs

    gcutm = 2*pi*maximum( pw.Ns )  
    #gcutm = pw.ecutrho
    #gcutm = maximum(pw.gvec.G2)/2
    @printf("gcutm, ecutrho, max G2 = %18.10f, %18.10f, %18.10f\n", gcutm, pw.ecutrho, maximum(pw.gvec.G2))

    alat = 1.0
    gamma_only = false
    gstart = 2

    charge = 0.0
    for ia = 1:Natoms
        isp = atm2species[ia]
        charge = charge + Zv[isp]
    end
    
    alpha = 2.9

    do_loop_alpha = true
    while(do_loop_alpha)
        alpha = alpha - 0.1
        #
        # choose alpha in order to have convergence in the sum over G
        # upperbound is a safe upper bound for the error in the sum over G
        if alpha <= 0.0
            @printf("ERROR in calculating Ewald energy:\n")
            @printf("optimal alpha not found\n")
            exit()
        end
        #
        # beware of unit of gcutm
        upperbound = 2.0*charge^2 * sqrt(2.0*alpha/2.0/pi) * erfc(sqrt(gcutm/4.0/alpha))  
        @printf("alpha, upperbound = %f %18.10e\n", alpha, upperbound)
        if upperbound <= 1.e-7
            @printf("do_loop_alpha is false, alpha = %18.10f\n", alpha)
            do_loop_alpha = false
        end
    end
    @printf("alpha = %18.10f\n", alpha)
    #
    # G-space sum here.
    # Determine if this processor contains G=0 and set the constant term
    #
    if gstart==2
        ewaldg = -charge^2 / alpha / 4.0
    else
        ewaldg = 0.0
    end

    # gamma_only should be .FALSE. for our case
    if gamma_only
        fact = 2.0
    else
        fact = 1.0
    end

    for ig = gstart:Ng
        rhon = 0.0 + im*0.0
        for isp = 1:Nspecies
            rhon = rhon + Zv[isp]*conj(Sf[ig,isp])
        end
        ewaldg = ewaldg + fact*abs(rhon)^2 * exp( -G2[ig]/alpha/4.0 ) / G2[ig]
    end
    ewaldg = 2.0 * 2*pi / omega * ewaldg
    #
    #  Here add the other constant term
    #
    if gstart==2
        for ia = 1:Natoms
            isp = atm2species[ia]
            ewaldg = ewaldg - Zv[isp]^2 * sqrt(8.0/2.0/pi*alpha)
        end
    end
    #
    # R-space sum here (only for the processor that contains G=0)
    #
    mxr = 50
    ewaldr = 0.0
    if gstart == 2
        rmax = 4.0 / sqrt(alpha) / alat
        #
        # with this choice terms up to ZiZj*erfc(4) are counted (erfc(4)=2x10^-8
        #
        for ia = 1:Natoms
            for ja = 1:Natoms
                dtau[:] = atpos[:,ia] - atpos[:,ja]
                # generates nearest-neighbors shells
                nrm, r, r2 = rgen( dtau, rmax, mxr, at, bg )
                # and sum to the real space part
                for ir = 1:nrm
                    rr = sqrt( r2[ir] ) * alat
                    isp = atm2species[ia]
                    jsp = atm2species[ja]
                    ewaldr = ewaldr + Zv[isp] * Zv[jsp] * erfc(sqrt(alpha)*rr)/rr
                end
            end
        end
    end
    
    @printf("ewaldg, ewaldr = %18.10f %18.10f\n", ewaldg, ewaldr)
    E_nn = 0.5*(ewaldg + ewaldr)
    return E_nn

end


