
#
# Based on Ylm calculation used in PWSCF
#
function Ylm_real( l::Int64, m::Int64, R::Array{Float64,1} )

    SMALL = 1.0e-9

    Rmod = sqrt( R[1]^2 + R[2]^2 + R[3]^2 )
    if Rmod < SMALL
        cost = 0.0
    else
        cost = R[3]/Rmod
    end

    #
    # beware the arc tan, it is defined modulo pi
    #
    if R[1] > SMALL
        phi = atan( R[2]/R[1] )
    elseif R[1] < -SMALL
        phi = atan( R[2]/R[1] ) + pi
    else
        #phi = pi/2.0*sign(R[2])  # XXX If R[2] == 0 ???
        phi = if R[2] >= 0 pi/2 else -pi/2 end
    end
    sint = sqrt( max(0.0, 1.0 - cost^2) )


    ylm = 0.0

    if l == 0
        ylm = 0.5*sqrt(1.0/pi)
        return ylm

    elseif l == 1
        # py
        if m == -1
            ylm = 0.5*sqrt(3.0/pi)*sint*sin(phi)
            return ylm
        # pz
        elseif m == 0
            ylm = 0.5*sqrt(3.0/pi)*cost
            return ylm
        # px
        elseif m == 1
            ylm = 0.5*sqrt(3.0/pi)*sint*cos(phi)
            return ylm
        end

    elseif l == 2
        # dxy
        if m == -2
            ylm = sqrt(15.0/16.0/pi) * sint^2 * sin(2.0*phi)
            return ylm
        # dyz
        elseif m == -1
            ylm = sqrt(15.0/4.0/pi)*cost*sint*sin(phi)
            return ylm
        # dz2
        elseif m == 0
            ylm = 0.25*sqrt(5.0/pi)*( 3.0*cost^2 - 1.0 )
            return ylm
        # dxz
        elseif m == 1
            ylm = sqrt(15.0/4.0/pi)*cost*sint*cos(phi)
            return ylm
        # dx2-y2
        elseif m == 2
            ylm = sqrt(15.0/16.0/pi) * sint^2 * cos(2.0*phi)
            return ylm
        end

    elseif l == 3

        if m == -3
            ylm = 0.25*sqrt(35.0/2.0/pi) * sint^3 * sin(3.0*phi)
            return ylm

        elseif m == -2
            ylm = 0.25*sqrt(105.0/pi)* sint^2 *cost * sin(2.0*phi)
            return ylm

        elseif m == -1
            ylm = 0.25*sqrt(21.0/2.0/pi)*sint*( 5.0*cost^2 - 1.0 )*sin(phi)
            return ylm

        elseif m == 0
            ylm = 0.25*sqrt(7.0/pi)*( 5.0*cost^3 - 3.0*cost )
            return ylm

        elseif m == 1
            ylm = 0.25*sqrt(21.0/2.0/pi)*sint*( 5.0*cost^2 - 1.0 )*cos(phi)
            return ylm

        elseif m == 2
            ylm = 0.25*sqrt(105.0/pi) * sint^2 * cost * cos(2.0*phi)
            return ylm

        elseif m == 3
            ylm = 0.25*sqrt(35.0/2.0/pi) * sint^3 * cos(3.0*phi)
            return ylm
        end

    end

end
