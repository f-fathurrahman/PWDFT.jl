using Printf
using FFTW

#import Plots
#Plots.gr()

#Plots.default( markersize=4 )
#Plots.default( markershape=:circle )
#Plots.default( framestyle=:box )
#Plots.default( leg=false )

# a periodic function
function funcx( cx, L, x, alpha )
    return exp( alpha*cos(2.0*pi/L *(x - cx) ) )
end


function test_main()
    Npoints = 25
    L = 10.0
    x = range(0.0, stop=L, length=Npoints)
    fx = zeros(Npoints)
    cx = 5.0
    alpha = 1.0
    for i = 1:Npoints
        fx[i] = funcx( cx, L, x[i], alpha)
    end
    
    #Plots.plot(x, fx)
    #Plots.savefig("TEMP_funcx_1d.pdf")

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

    for ip = 1:Npoints
        s = 0.0 + im*0.0
        for ig = 1:Ng
            ipp = idx_g2r[ig]
            s = s + fg[ipp]*exp(im*G[ig]*x[ip])/sqrt(L)
        end
        @printf("%4d %18.10f %18.10f %18.10f\n", ip, fx[ip], real(s), imag(s))
    end
end

function mm_to_nn( mm::Int64,S::Int64 )
    if mm > S/2
        return mm - S
    else
        return mm
    end
end

test_main()