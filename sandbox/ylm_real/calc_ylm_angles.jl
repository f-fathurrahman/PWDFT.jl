
#
# Based on Ylm calculation used in PWSCF
#
function calc_ylm_angles( R::Array{Float64,1} )

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

    @printf("phi  = %18.10f\n", phi)
    @printf("cost = %18.10f\n", cost)
    @printf("sint = %18.10f\n", sint)

    return

end
