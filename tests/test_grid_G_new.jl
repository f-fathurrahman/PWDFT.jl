using PWDFT

function mm_to_nn(mm::Int64,S::Int64)
    if mm > S/2
        return mm - S
    else
        return mm
    end
end


function calc_Ng( Ns, RecVecs, ecutrho )

    ig = 0
    Ngs = 0
    G = zeros(Float64,3)

    for k in 0:Ns[3]-1
    for j in 0:Ns[2]-1
    for i in 0:Ns[1]-1
        ig = ig + 1
        gi = mm_to_nn( i, Ns[1] )
        gj = mm_to_nn( j, Ns[2] )
        gk = mm_to_nn( k, Ns[3] )
        G[1] = RecVecs[1,1]*gi + RecVecs[2,1]*gj + RecVecs[3,1]*gk
        G[2] = RecVecs[1,2]*gi + RecVecs[2,2]*gj + RecVecs[3,2]*gk
        G[3] = RecVecs[1,3]*gi + RecVecs[2,3]*gj + RecVecs[3,3]*gk
        G2 = G[1]^2 + G[2]^2 + G[3]^2
        if 0.5*G2 < ecutrho
            Ngs = Ngs + 1
        end
    end
    end
    end

    return Ngs

end

struct GVectorsNew
    Ng::Int64
    G::Array{Float64,2}
    G2::Array{Float64,1}
    idx_g2r::Array{Int64,1}
end


function init_grid_G_new( Ns, RecVecs, ecutrho )

    Ng = calc_Ng( Ns, RecVecs, ecutrho )

    G_temp = zeros(Float64,3)

    G  = Array{Float64}(3,Ng)
    G2 = Array{Float64}(Ng)
    idx_g2r = Array{Int64}(Ng)

    ig = 0
    ip = 0
    for k in 0:Ns[3]-1
    for j in 0:Ns[2]-1
    for i in 0:Ns[1]-1
        ip = ip + 1
        gi = mm_to_nn( i, Ns[1] )
        gj = mm_to_nn( j, Ns[2] )
        gk = mm_to_nn( k, Ns[3] )
        G_temp[1] = RecVecs[1,1]*gi + RecVecs[2,1]*gj + RecVecs[3,1]*gk
        G_temp[2] = RecVecs[1,2]*gi + RecVecs[2,2]*gj + RecVecs[3,2]*gk
        G_temp[3] = RecVecs[1,3]*gi + RecVecs[2,3]*gj + RecVecs[3,3]*gk
        G2_temp = G_temp[1]^2 + G_temp[2]^2 + G_temp[3]^2
        if 0.5*G2_temp < ecutrho
            ig = ig + 1
            G[:,ig] = G_temp[:]
            G2[ig] = G2_temp
            idx_g2r[ig] = ip
        end
    end
    end
    end

    return GVectorsNew( Ng, G, G2, idx_g2r )

end

import PWDFT.Poisson_solve
function Poisson_solve( pw::PWGrid, gvec::GVectorsNew, rhoR )
    #
    G2 = gvec.G2
    Ng = gvec.Ng
    idx_g2r = gvec.idx_g2r

    Ns = pw.Ns
    Npoints = prod(Ns)
    #
    ctmp = 4.0*pi*R_to_G( pw, rhoR )
    ctmp[1] = 0.0
    for ig = 2:Ng
        ip = idx_g2r[ig]
        ctmp[ip] = ctmp[ip]/G2[ig]
    end
    return ctmp
end

function test_main()
    LatVecs = 16.0*diagm(ones(3))
    ecutwfc_Ry = 40.0
    pw = PWGrid( ecutwfc_Ry*0.5, LatVecs )
    println(pw)

    Ns = pw.Ns
    RecVecs = pw.RecVecs
    ecutrho = pw.ecutrho

    gvec = init_grid_G_new( Ns, RecVecs, ecutrho )
    println("Ng = ", gvec.Ng)
end

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


function test_Poisson( ecutwfc_Ry::Float64 )
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
    phi = real( G_to_R(pw, phiG) )
    Ehartree = 0.5*dot( phi, rho ) * Ω/Npoints
    Uanal = ( (1/σ1 + 1/σ2)/2 - sqrt(2) / sqrt( σ1^2 + σ2^2 ) ) / sqrt(pi)
    @printf("Num, ana, diff = %18.10f %18.10f %18.10e\n", Ehartree, Uanal, abs(Ehartree-Uanal))
    #
    # Using GvectorsNew
    #
    Ns = pw.Ns
    RecVecs = pw.RecVecs
    ecutrho = pw.ecutrho
    gvec_new = init_grid_G_new( Ns, RecVecs, ecutrho )
    println("Ng = ", gvec_new.Ng)    
    #
    phiG_new = Poisson_solve( pw, gvec_new, rho )
    phi_new = real( G_to_R(pw, phiG_new) )
    Ehartree = 0.5*dot( phi_new, rho ) * Ω/Npoints
    Uanal = ( (1/σ1 + 1/σ2)/2 - sqrt(2) / sqrt( σ1^2 + σ2^2 ) ) / sqrt(pi)
    @printf("Num, ana, diff = %18.10f %18.10f %18.10e\n", Ehartree, Uanal, abs(Ehartree-Uanal))    
end

#test_main()

test_Poisson(30.0)