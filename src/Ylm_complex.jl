function Ylm_complex( l::Int64, m::Int64, R::Array{Float64,1} )

    const SMALL = 1.0e-9

    x = R[1]
    y = R[2]
    z = R[3]

    r = sqrt( x^2 + y^2 + z^2 )
    if r < SMALL
        return 0.0 + im*0.0
    end

    ylm = 0.0 + im*0.0

    if l == 0
        ylm = 0.5*sqrt(1.0/pi) + im*0.0
        return ylm

    elseif l == 1
        # py
        if m == -1
            ylm = 0.5*sqrt(3.0/2.0/pi)*( x - im*y )/r
            return ylm
        # pz
        elseif m == 0
            ylm = 0.5*sqrt(3.0/pi)*z/r + im*0.0
            return ylm
        # px
        elseif m == 1
            ylm = -0.5*sqrt(3.0/2.0/pi)*( x + im*y )/r
            return ylm
        end

    elseif l == 2
        # dxy
        if m == -2
            ylm = 0.25*sqrt(15.0/2.0/pi) * ( x - im*y )^2 / r^2
            return ylm
        # dyz
        elseif m == -1
            ylm = 0.5*sqrt(15.0/2.0/pi) * ( x - im*y ) * z / r^2
            return ylm
        # dz2
        elseif m == 0
            ylm = 0.25*sqrt(5.0/pi) * (2.0*z^2 - x^2 - y^2) / r^2 + im*0.0
            return ylm
        # dxz
        elseif m == 1
            ylm = -0.5*sqrt(15.0/2.0/pi) * ( x + im*y ) * z / r^2
            return ylm
        # dx2-y2
        elseif m == 2
            ylm = 0.25*sqrt(15.0/2.0/pi) * ( x + im*y )^2 / r^2
            return ylm
        end

    elseif l == 3

        if m == -3
            ylm = 0.125*sqrt(35.0/pi) * ( x - im*y )^3 / r^3
            return ylm

        elseif m == -2
            ylm = 0.25*sqrt(105.0/2.0/pi) * ( x - im*y )^2 * z / r^3
            return ylm

        elseif m == -1
            ylm = 0.125 * sqrt(21.0/pi) * ( x - im*y ) * ( 4.0*z^2 - x^2 - y^2 ) / r^3
            return ylm

        elseif m == 0
            ylm = 0.25*sqrt(7.0/pi) * ( 2.0*z^2 - 3.0*x^2 - 3*y^2 ) / r^3 + im*0.0
            return ylm

        elseif m == 1
            ylm = -0.125 * sqrt(21.0/pi) * ( x + im*y ) * ( 4.0*z^2 - x^2 - y^2 ) / r^3
            return ylm

        elseif m == 2
            ylm = 0.25*sqrt(105.0/2.0/pi) * ( x + im*y )^2 * z / r^3
            return ylm

        elseif m == 3
            ylm = -0.125*sqrt(35.0/pi) * ( x + im*y )^3 / r^3
            return ylm
        end

    end


end
