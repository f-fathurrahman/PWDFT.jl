#
# Based on Ylm calculation used in PWSCF
#
function Ylm_real_qe!(
    lmax::Int64,
    R::Matrix{Float64},
    Ylm::Matrix{Float64}
)

    @assert size(R,1) == 3

    N = size(R,2)
    @assert size(Ylm,1) == N

    lmmax = (lmax + 1)^2
    @assert size(Ylm,2) == lmmax

    for i in 1:N
        @views Ylm_real_qe!(lmax, R[:,i], Ylm[i,:])
    end
    return
end


function Ylm_real_qe!(
    lmax::Int64,
    R::AbstractVector{Float64},
    Ylm::AbstractVector{Float64}
)

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
        if R[2] >= 0
            phi = pi/2
        else
            phi = -pi/2
        end
    end
    sint = sqrt( max(0.0, 1.0 - cost^2) )

    Q = OffsetArray( zeros(lmax+1,lmax+1), 0:lmax, 0:lmax )

    #  Q(:,l,m) are defined as sqrt((l-m)!/(l+m)!) * P(:,l,m) where
    #  P(:,l,m) are the Legendre Polynomials (0 <= m <= l)
    lm = 0
    for l in 0:lmax
        c = sqrt((2*l+1)/4π)
        if l == 0
            Q[0,0] = 1.0
        elseif l == 1
            Q[1,0] =  cost
            Q[1,1] = -sint/sqrt(2)
        else
            # recursion on l for Q(:,l,m)
            for m in 0:(l-2)
                Q[l,m] = cost*(2*l-1)/sqrt(l*l - m*m) * Q[l-1,m] -
                         sqrt((l-1)*(l-1) - m*m) / sqrt(l*l - m*m) * Q[l-2,m]
            end
            #
            Q[l,l-1] = cost * sqrt(2*l-1) * Q[l-1,l-1]
            Q[l,l] = -sqrt(2*l - 1)/sqrt(2*l) * sint * Q[l-1,l-1]
        end
        # Y_lm, m = 0
        lm = lm + 1
        Ylm[lm] = c*Q[l,0]
        for m in 1:l
            # Y_lm, m > 0
            lm = lm + 1
            Ylm[lm] = c * sqrt(2) * Q[l,m] * cos(m*phi)
            # Y_lm, m < 0
            lm = lm + 1
            Ylm[lm] = c * sqrt(2) * Q[l,m] * sin(m*phi)
        end
    end
    return
end

