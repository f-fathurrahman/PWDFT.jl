using Printf
using FFTW

import Plots
Plots.gr()

Plots.default( markersize=4 )
Plots.default( markershape=:circle )
Plots.default( framestyle=:box )
Plots.default( leg=false )

# a periodic function
function funcx( cx, L, x, alpha )
    return exp( alpha*cos(2.0*pi/L *(x - cx) ) )
end

function gauss1d( alpha::Float64, sigma::Float64, r::Float64 )
    r2 = r[1]^2
    c1 = 2*sigma^2
    cc1 = sqrt(2*pi*sigma^2)
    return alpha*exp(-r2/c1)/cc1
end

function test_main()
    Npoints = 25
    L = 10.0
    x = range(0.0, stop=L, length=Npoints)
    fx = zeros(Npoints)
    cx = 5.0
    alpha = 1.0
    sigma = 1.0
    for ip = 1:Npoints
        #fx[ip] = funcx( cx, L, x[ip], alpha)
        fx[ip] = gauss1d(alpha, sigma, x[ip]-cx)
    end

    fg = fft( fx )

    Ng = Npoints
    G = zeros(Ng)
    ig = 0
    idx_g2r = zeros(Int64,Ng)
    for ip = 1:Npoints
        ii = mm_to_nn(ip,Npoints)
        ig = ig + 1
        G[ig] = ii*2*pi/L
        idx_g2r[ig] = ip
        @printf("ip = %4d, ig = %4d, ii=%4d\n", ip, ig, ii)
    end

    println()

    f_interp = zeros(ComplexF64,Ng)
    for ip = 1:Npoints
        s = 0.0 + im*0.0
        for ig = 1:Ng
            ipp = idx_g2r[ig]
            s = s + fg[ipp]*exp(im*G[ig]*x[ip])/L
        end
        f_interp[ip] = s
        @printf("%4d %18.10f %18.10f %18.10f\n", ip, fx[ip], real(s), imag(s))
    end

    Plots.plot(x, abs.(real(f_interp)), leg=true)
    Plots.plot!(x, imag(f_interp))
    Plots.plot!(x, abs.(f_interp))
    Plots.plot!(x, fx)
    Plots.savefig("TEMP_funcx_1d_interp.pdf")
end

function mm_to_nn( mm::Int64,S::Int64 )
    if mm > S/2
        return mm - S
    else
        return mm
    end
end

test_main()