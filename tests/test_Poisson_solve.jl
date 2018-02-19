using PWDFT
using PWDFT.PW

include("Poisson_solve.jl")

function gen_dr( r, center )
    Npoints = size(r)[2]
    dr = Array{Float64}(Npoints)
    #
    for ip=1:Npoints
        dx2 = ( r[1,ip] - center[1] )^2
        dy2 = ( r[2,ip] - center[2] )^2
        dz2 = ( r[3,ip] - center[3] )^2
        dr[ip] = sqrt( dx2 + dy2 + dz2 )
    end
    return dr
end

function gen_rho( dr, σ1, σ2 )
    Npoints = size(dr)[1]
    rho = Array{Float64}(Npoints)
    c1 = 2*σ1^2
    c2 = 2*σ2^2
    cc1 = sqrt(2*pi*σ1^2)^3
    cc2 = sqrt(2*pi*σ2^2)^3
    for ip=1:Npoints
        g1 = exp(-dr[ip]^2/c1)/cc1
        g2 = exp(-dr[ip]^2/c2)/cc2
        rho[ip] = g2 - g1
    end
    return rho
end


function test_main( ecutwfc_Ry::Float64 )
    #
    LatVecs = 16.0*diagm( ones(3) )
    #
    pw = PWGrid( ecutwfc_Ry*0.5, LatVecs )
    println(pw)
    #
    Npoints = prod(pw.Ns)
    Ω = pw.Ω
    r = pw.r
    Ns = pw.Ns
    #
    # Generate array of distances
    #
    center = sum(LatVecs,2)/2
    dr = gen_dr( pw.r, center )
    #
    # Generate charge density
    #
    const σ1 = 0.75
    const σ2 = 0.50
    rho = gen_rho( dr, σ1, σ2 )
    #
    # Solve Poisson equation and calculate Hartree energy
    #
    phiG = Poisson_solve( pw, rho )
    phi = real( c_G_to_R(Ns, phiG) )
    Ehartree = 0.5*dot( phi, rho ) * Ω/Npoints
    #
    Uanal = ( (1/σ1 + 1/σ2)/2 - sqrt(2) / sqrt( σ1^2 + σ2^2 ) ) / sqrt(pi)
    @printf("Num, ana, diff = %18.10f %18.10f %18.10e\n", Ehartree, Uanal, abs(Ehartree-Uanal))
end

@code_native test_main(1.0)
@time test_main(30.0)
