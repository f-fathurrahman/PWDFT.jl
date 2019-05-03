
#
# Based on Ylm calculation used in PWSCF
#
function Ylm_real_v2( l::Int64, m::Int64, R::Array{Float64,1} )

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


    fpi = 4.0*pi
    c = sqrt( (2*l+1)/fpi )

    if m == 0
        return c * eval_Q_lm( l, 0, cost, sint)

    elseif m > 0
        return c*sqrt(2.0) * eval_Q_lm( l, m, cost, sint) * cos( m*phi )

    elseif m < 0
        return c*sqrt(2.0) * eval_Q_lm( l, -m, cost, sint) * sin( -m*phi )

    end

end


function eval_Q_lm(l, m, cost, sint)

    prefactor = sqrt( factorial(l-m) / factorial(l+m) )
    
    Q = 0.0

    if l == 0
        return Q = 1.0
    
    elseif l == 1
        if m == 0
            Q = cost
        elseif m == 1
            Q = -sint
        else
            error("Should not reached here")
        end
    
    elseif l == 2
        if m == 0
            Q = 0.5*(3*cost^2 - 1)
        elseif m == 1
            Q = -3.0*cost*sint
        elseif m == 2
            Q = 3*sint^2
        else
            error("Should not reached here")
        end

    elseif l == 3
        if m == 0
            Q = 0.5*(5*cost^3 - 3*cost)
        elseif m == 1
            Q = -1.5*(5*cost^2 - 1)*sint
        elseif m == 2
            Q = 15*cost*sint^2
        elseif m == 3
            Q = -15*sint^3
        else
            error("Should not reached here")
        end

    else

        error("Should not reached here")

    end

    return prefactor*Q
end