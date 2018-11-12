using Printf
using PWDFT

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

# max: alpha
function gauss3d( alpha::Float64, sigma::Float64, r::Tuple{Float64,Float64,Float64} )
    r2 = r[1]^2 + r[2]^2 + r[3]^2
    c1 = 2*sigma^2
    cc1 = sqrt(2*pi*sigma^2)^3
    return alpha*exp(-r2/c1)/cc1
end

function plot_funcx()
    L = 10.0
    Lplot = 2*L
    Nx = 100
    cx1 = L/3
    cx2 = 2*L/3
    alpha1 = 1.0
    alpha2 = 1.0
    x = range( 0.0, stop=Lplot, length=100 )
    y = zeros(Nx)
    for i = 1:Nx
        y[i] = funcx( cx1, L, x[i], alpha1 ) + funcx( cx2, L, x[i], alpha2 )
    end
    Plots.plot(x, y)
    Plots.savefig("TEMP_funcx.pdf")
end

function test_main()
    L = 10.0
    ecutwfc = 25.0
    # choose a rather large cell
    LatVecs = gen_lattice_sc(L)
    
    pw = PWGrid( ecutwfc, LatVecs )
    println(pw)

    Npoints = prod(pw.Ns)
    dVol = pw.CellVolume/Npoints
    fr = zeros(Npoints)
    c = (5.0, 5.0, 5.0)
    alpha = 1.0
    sigma = 0.3
    L = 10.0
    for ip = 1:Npoints
        x = pw.r[1,ip]
        y = pw.r[2,ip]
        z = pw.r[3,ip]
        #fr[ip] = funcx(c[1], L, x, alpha)*
        #         funcx(c[2], L, y, alpha)*
        #         funcx(c[3], L, z, alpha)
        r = (x-c[1], y-c[2], z-c[3])
        fr[ip] = gauss3d( alpha, sigma, r )
    end
    integ_fr = sum(fr)*dVol
    @printf("integ_fr  = %18.10f\n", integ_fr)

    fg = R_to_G( pw, fr )  # return all

    fgr = real(G_to_R( pw, fg ))
    integ_fgr = sum(fgr)*dVol
    @printf("integ_fgr = %18.10f\n", integ_fgr)

    # scale with structure factor
    #=G = pw.gvec.G
    Ng = pw.gvec.Ng
    for ig = 1:Ng
        ip = pw.gvec.idx_g2r[ig]
        GX = c[1]*G[1,ig] + c[2]*G[2,ig] + c[3]*G[3,ig]
        Sf = cos(GX) - im*sin(GX)
        fg[ip] = fg[ip]*Sf
    end
    =#

    r_eval = (5.0, 5.0, 5.0)
    rr = ( r_eval[1]-c[1], r_eval[2]-c[2], r_eval[3]-c[3] )
    f_interp = pw_interp(pw, fg, r_eval)
    @printf("fr       = %18.10f\n", gauss3d(alpha, sigma, rr) )
    @printf("f_interp = %18.10f\n", real(f_interp)/Npoints)

    println("ratio = ", f_interp/Npoints/gauss3d(alpha, sigma, rr))
    

    #=
    f_interp = zeros(Npoints)
    for ip = 1:Npoints
        rr = pw.r[:,ip]
        f_interp[ip] = real( pw_interp, rr )
    end
    =#
end


function pw_interp( pw::PWGrid, fg::Array{ComplexF64,1}, r::Tuple{Float64,Float64,Float64} )
    s = 0.0 + im*0.0
    for ig = 1:pw.gvec.Ng
        ip = pw.gvec.idx_g2r[ig]
        G = @view pw.gvec.G[:,ig]
        Gr = G[1]*r[1] + G[2]*r[2] + G[3]*r[3]
        s = s + fg[ip]*exp( im*Gr )
    end
    return s
end

#plot_funcx()
test_main()