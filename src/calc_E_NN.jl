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

"""
Calculates Ewald energy with both G- and R-space terms.
Determines optimal alpha. Should hopefully work for any structure.
"""
function calc_E_NN( pw::PWGrid, atoms::Atoms, Zv::Array{Float64} )

    Sf = calc_strfact( atoms, pw )

    atpos = atoms.positions
    Natoms = atoms.Natoms
    Nspecies = atoms.Nspecies
    atm2species = atoms.atm2species

    Ns = pw.Ns
    LL = pw.LatVecs
    Ng = pw.gvec.Ng
    G2 = pw.gvec.G2

    ewaldg = 0.0
    ewaldsr = 0.0

    dtau = zeros(Float64,3)

  REAL(DP) :: charge, ewaldg, ewaldr, dtau (3), alpha, &
       r (3, mxr), r2 (mxr), rmax, rr, upperbound, fact

    at = pw.LatVecs'  # assuming transpose

    omega = pw.Ω

    bg = pw.RecVecs

    gcutm = 2*pi*maxval( pw.Ns )  

    alat = 1.d0
    gamma_only = .FALSE.
    gstart = 2

    charge = 0.d0
    for ia = 1:Natoms
        isp = atm2species[ia]
        charge = charge + Zv[isp]
    end
    
    alpha = 2.9

    do_loop_alpha = true
    while(do_loop_alpha)
        alpha = alpha - 0.1d0
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
        upperbound = 2.0*charge^2 * sqrt(2.0*alpha/(2*pi)) * erfc(sqrt(gcutm/4.0/alpha))  
        if upperbound <= 1.e-7
            do_loop_alpha = false
        end
    end
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

    if ng = gstart:ngm
        rhon = 0.0 + im*0.0
        for isp = 1:Nspecies
            rhon = rhon + Zv[nt]*conj(strf[ng,isp])
        end
        ewaldg = ewaldg + fact*abs(rhon)^2 * exp( -gg[ng]/alpha/4.d0 )/ gg[ng]
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
    ewaldr = 0.d0
    if gstart==2
        rmax = 4.d0 / sqrt(alpha) / alat
        #
        # with this choice terms up to ZiZj*erfc(4) are counted (erfc(4)=2x10^-8
        #
        for na = 1:Natoms
            for nb = 1:Natoms
                dtau[:] = tau[:,na] - tau[:,nb]
                # generates nearest-neighbors shells
                CALL rgen( dtau, rmax, mxr, at, bg, r, r2, nrm )
                # and sum to the real space part
                for nr = 1:nrm
                    rr = sqrt( r2[nr] ) * alat
                    ewaldr = ewaldr + Zv[ityp[na]] * Zv[ityp[nb]] * erfc(sqrt(alpha)*rr)/rr
                end
            end
        end
    end
  
    E_nn = 0.5d0*(ewaldg + ewaldr)

end


function rgen( dtau, rmax, mxr, at, bg )
    #
    #   generates neighbours shells (cartesian, in units of lattice parameter)
    #   with length < rmax,and returns them in order of increasing length:
    #      r(:) = i*a1(:) + j*a2(:) + k*a3(:) - dtau(:),   r2 = r^2
    #   where a1, a2, a3 are primitive lattice vectors. Other input variables:
    #     mxr = maximum number of vectors
    #     at  = lattice vectors ( a1=at(:,1), a2=at(:,2), a3=at(:,3) )
    #     bg  = reciprocal lattice vectors ( b1=bg(:,1), b2=bg(:,2), b3=bg(:,3) )
    #   Other output variables:
    #     nrm = the number of vectors with r^2 < rmax^2
    #

    nrm = 0
    r = zeros(Float64,3,mxr)
    r2 = zeros(Float64,mxr)

  INTEGER, ALLOCATABLE :: irr (:)
  INTEGER ::  nm1, nm2, nm3, i, j, k, ipol, ir, indsw, iswap
  real(DP) :: tt, swap
  real(DP), EXTERNAL :: dnrm2

  ds = zeros(Float64,3)
  dtau0 = zeros(Float64,3)
  t = zeros(Float64,3)
  !
  !
  nrm = 0
  IF (rmax==0.d0) RETURN

  ! bring dtau into the unit cell centered on the origin - prevents trouble
  ! if atomic positions are not centered around the origin but displaced
  ! far away (remember that translational invariance allows this!)
  !
  ds(:) = matmul( dtau(:), bg(:,:) )
  ds(:) = ds(:) - anint(ds(:))
  dtau0(:) = matmul( at(:,:), ds(:) )
  !
  ALLOCATE (irr( mxr))
  !
  ! these are estimates of the maximum values of needed integer indices
  !
  nm1 = int (dnrm2 (3, bg (1, 1), 1) * rmax) + 2
  nm2 = int (dnrm2 (3, bg (1, 2), 1) * rmax) + 2
  nm3 = int (dnrm2 (3, bg (1, 3), 1) * rmax) + 2
  !
  DO i = -nm1, nm1
     DO j = -nm2, nm2
        DO k = -nm3, nm3
           tt = 0.d0
           DO ipol = 1, 3
              t (ipol) = i*at (ipol, 1) + j*at (ipol, 2) + k*at (ipol, 3) &
                       - dtau0(ipol)
              tt = tt + t (ipol) * t (ipol)
           ENDDO
           IF (tt<=rmax**2.and.abs (tt) >1.d-10) THEN
              nrm = nrm + 1
              IF (nrm>mxr) then 
                WRITE(*,*) 'ERROR in rgen: too many r-vectors', nrm
                STOP 
              ENDIF 
              DO ipol = 1, 3
                 r (ipol, nrm) = t (ipol)
              ENDDO
              r2 (nrm) = tt
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  !
  !   reorder the vectors in order of increasing magnitude
  !
  !   initialize the index inside sorting routine
  !
  irr (1) = 0
  IF (nrm>1) CALL hpsort (nrm, r2, irr)
  DO ir = 1, nrm - 1
20   indsw = irr (ir)
     IF (indsw/=ir) THEN
        DO ipol = 1, 3
           swap = r (ipol, indsw)
           r (ipol, indsw) = r (ipol, irr (indsw) )
           r (ipol, irr (indsw) ) = swap
        ENDDO
        iswap = irr (ir)
        irr (ir) = irr (indsw)
        irr (indsw) = iswap
        GOTO 20
     ENDIF

  ENDDO
  DEALLOCATE(irr)
  !
  RETURN
END SUBROUTINE rgen
