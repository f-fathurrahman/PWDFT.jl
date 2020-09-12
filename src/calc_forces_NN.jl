using SpecialFunctions: erfc

function calc_forces_NN( atoms::Atoms )
    F_NN = zeros(3, atoms.Natoms)
    calc_forces_NN!( atoms, F_NN )
    return F_NN
end

function calc_forces_NN!( atoms::Atoms, F_NN::Array{Float64,2} )
    return calc_forces_NN!( atoms.LatVecs, atoms, atoms.Zvals, F_NN )
end

function _Herfc(x)
    return -2*exp(-x^2)/sqrt(pi) - erfc(x)/x
end

function calc_forces_NN!(
    LatVecs::Array{Float64,2},
    atoms::Atoms,
    Zvals::Array{Float64,1},
    F_NN::Array{Float64,2}
)

    t1 = LatVecs[:,1]
    t2 = LatVecs[:,2]
    t3 = LatVecs[:,3]
  
    Ω = abs(det(LatVecs))

    RecVecs = 2*pi*inv(LatVecs')
    g1 = RecVecs[:,1]
    g2 = RecVecs[:,2]
    g3 = RecVecs[:,3]

    t1m = sqrt(dot(t1,t1))
    t2m = sqrt(dot(t2,t2))
    t3m = sqrt(dot(t3,t3))

    g1m = sqrt(dot(g1,g1))
    g2m = sqrt(dot(g2,g2))
    g3m = sqrt(dot(g3,g3))

    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    # Atomic positions
    tau = atoms.positions

    # Parameters
    gcut = 20.0
    ebsl = 1e-10

    glast2 = gcut*gcut
    gexp = -log(ebsl)    
    η = sqrt(glast2/gexp)/2

    tmax = sqrt(0.5*gexp)/η

    mmm1 = round(Int64, tmax/t1m + 1.5)
    mmm2 = round(Int64, tmax/t2m + 1.5)
    mmm3 = round(Int64, tmax/t3m + 1.5)

    #println("R-space: tmax = ", tmax)
    #println("R-space: mmm1 = ", mmm1)
    #println("R-space: mmm2 = ", mmm2)
    #println("R-space: mmm3 = ", mmm3)

    dtau = zeros(Float64,3)
    G = zeros(Float64,3)
    T = zeros(Float64,3)

    F_NN_R = zeros(3,Natoms)
    F_NN_G = zeros(3,Natoms)    

    D = zeros(3)

    for ia in 1:Natoms, ja in ia+1:Natoms

        #@printf("ia = %d ja = %d\n", ia, ja)

        dtau[1] = tau[1,ia] - tau[1,ja]
        dtau[2] = tau[2,ia] - tau[2,ja]
        dtau[3] = tau[3,ia] - tau[3,ja]
        
        isp = atm2species[ia]
        jsp = atm2species[ja]
        ZiZj = Zvals[isp]*Zvals[jsp]

        for i in -mmm1:mmm1, j in -mmm2:mmm2, k in -mmm3:mmm3
            if (abs(i) + abs(j) + abs(k)) == 0
                continue
            end
            T[1] = i*t1[1] + j*t2[1] + k*t3[1]
            T[2] = i*t1[2] + j*t2[2] + k*t3[2]
            T[3] = i*t1[3] + j*t2[3] + k*t3[3]
            D[1] = dtau[1] - T[1]
            D[2] = dtau[2] - T[2]
            D[3] = dtau[3] - T[3]
            D2 = D[1]^2 + D[2]^2 + D[3]^2
            Dmag = sqrt(D2)
            xx = _Herfc(η*Dmag)
            #
            df1 = η*xx*D[1]/D2*ZiZj
            df2 = η*xx*D[2]/D2*ZiZj
            df3 = η*xx*D[3]/D2*ZiZj
            #
            F_NN_R[1,ia] = F_NN_R[1,ia] - df1
            F_NN_R[2,ia] = F_NN_R[2,ia] - df2
            F_NN_R[3,ia] = F_NN_R[3,ia] - df3
            #
            F_NN_R[1,ja] = F_NN_R[1,ja] + df1
            F_NN_R[2,ja] = F_NN_R[2,ja] + df2
            F_NN_R[3,ja] = F_NN_R[3,ja] + df3
        end
    end

    mmm1 = round(Int64, gcut/g1m + 1.5)
    mmm2 = round(Int64, gcut/g2m + 1.5)
    mmm3 = round(Int64, gcut/g3m + 1.5)

    println("G-space: tmax = ", tmax)
    println("G-space: mmm1 = ", mmm1)
    println("G-space: mmm2 = ", mmm2)
    println("G-space: mmm3 = ", mmm3)

    for ia in 1:Natoms, ja in ia+1:Natoms
        
        #@printf("ia = %d ja = %d\n", ia, ja)

        isp = atm2species[ia]
        jsp = atm2species[ja]
        ZiZj = Zvals[isp]*Zvals[jsp]

        dtau[1] = tau[1,ia] - tau[1,ja]
        dtau[2] = tau[2,ia] - tau[2,ja]
        dtau[3] = tau[3,ia] - tau[3,ja]

        for i in -mmm1:mmm1, j in -mmm2:mmm2, k in -mmm3:mmm3
            if ( abs(i) + abs(j) + abs(k) ) == 0
                continue
            end
            G[1] = i*g1[1] + j*g2[1] + k*g3[1]
            G[2] = i*g1[2] + j*g2[2] + k*g3[2]
            G[3] = i*g1[3] + j*g2[3] + k*g3[3]        
            G2 = G[1]^2 + G[2]^2 + G[3]^2
            x = 4*pi/Ω * exp(-0.25*G2/η^2)/G2
            Gtau = G[1]*dtau[1] + G[2]*dtau[2] + G[3]*dtau[3]
            xsinGtau = x*sin(Gtau)
            df1 = xsinGtau*G[1]*ZiZj
            df2 = xsinGtau*G[2]*ZiZj
            df3 = xsinGtau*G[3]*ZiZj
            #
            F_NN_G[1,ia] = F_NN_G[1,ia] + df1
            F_NN_G[2,ia] = F_NN_G[2,ia] + df2
            F_NN_G[3,ia] = F_NN_G[3,ia] + df3
            #
            F_NN_G[1,ja] = F_NN_G[1,ja] - df1
            F_NN_G[2,ja] = F_NN_G[2,ja] - df2
            F_NN_G[3,ja] = F_NN_G[3,ja] - df3
        end
    end

    F_NN[:] = F_NN_G[:] + F_NN_R[:]

    return
end