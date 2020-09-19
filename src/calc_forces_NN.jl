function calc_forces_NN( pw::PWGrid, atoms::Atoms )
    F_NN = zeros(3, atoms.Natoms)
    calc_forces_NN!( pw, atoms, F_NN )
    return F_NN
end

function calc_forces_NN!( pw::PWGrid, atoms::Atoms, F_NN::Array{Float64,2} )
    return calc_forces_NN!( pw, atoms, atoms.Zvals, F_NN )
end

function _Herfc(x)
    return -2*exp(-x^2)/sqrt(pi) - erfc(x)/x
end

function calc_forces_NN!(
    pw::PWGrid,
    atoms::Atoms,
    Zvals::Array{Float64,1},
    F_NN::Array{Float64,2}
)

    LatVecs = pw.LatVecs  
    Ω = pw.CellVolume
    Natoms = atoms.Natoms
    atm2species = atoms.atm2species

    # Atomic positions
    tau = atoms.positions

    charge = 0.0
    for ia in 1:Natoms
        isp = atm2species[ia]
        charge = charge + Zvals[isp]
    end

    G2_max = pw.ecutrho
    α = 1.0
    upperbound = charge^2 * sqrt(α/π) * erfc(sqrt(G2_max/2.0/α))
    while upperbound > 1e-6
        α = α - 0.1
        if α <= 0.0
            error("Optimal α is not found")
        end
        upperbound = charge^2 * sqrt(α/π) * erfc(sqrt(G2_max/2.0/α))
    end

    Ng = pw.gvec.Ng
    G = pw.gvec.G
    G2 = pw.gvec.G2
    aux = zeros(ComplexF64,Ng)
    for ia in 1:Natoms
        isp = atm2species[ia]
        for ig in 2:Ng
            GX = tau[1,ia]*G[1,ig] + tau[2,ia]*G[2,ig] + tau[3,ia]*G[3,ig]
            Sf = cos(GX) + im*sin(GX) # conj
            aux[ig] = aux[ig] + Zvals[isp] * Sf
        end
    end

    # Treat 2d cutoff is skipped

    # skip G2=0
    for ig in 2:Ng
        aux[ig] = aux[ig] * exp(-G2[ig]/α/4.0) / G2[ig]
    end

    F_NN_G = zeros(3,Natoms)

    for ia in 1:Natoms
        isp = atm2species[ia]
        for ig in 2:Ng # from 1?
            GX = G[1,ig]*tau[1,ia] + G[2,ig]*tau[2,ia] + G[3,ig]*tau[3,ia]
            sumnb = cos(GX)*imag(aux[ig]) - sin(GX)*real(aux[ig])
            F_NN_G[1,ia] = F_NN_G[1,ia] + G[1,ig] * sumnb
            F_NN_G[2,ia] = F_NN_G[2,ia] + G[2,ig] * sumnb
            F_NN_G[3,ia] = F_NN_G[3,ia] + G[3,ig] * sumnb
        end
        fact = -4π * Zvals[isp]  / pw.CellVolume
        for i in 1:3
            F_NN_G[i,ia] = F_NN_G[i,ia]*fact
        end
    end

    # Real space sum

    dtau = zeros(Float64,3)
    T = zeros(Float64,3)
    F_NN_R = zeros(3,Natoms)
    D = zeros(3)

    η = sqrt(α)
    ebsl = 1e-10
    gexp = -log(ebsl)    
    tmax = sqrt(0.5*gexp)/η
    
    @views t1 = LatVecs[:,1]
    @views t2 = LatVecs[:,2]
    @views t3 = LatVecs[:,3]
    t1m = sqrt(dot(t1,t1))
    t2m = sqrt(dot(t2,t2))
    t3m = sqrt(dot(t3,t3))
    mmm1 = round(Int64, tmax/t1m + 1.5)
    mmm2 = round(Int64, tmax/t2m + 1.5)
    mmm3 = round(Int64, tmax/t3m + 1.5)

    for ia in 1:Natoms, ja in ia+1:Natoms

        dtau[1] = tau[1,ia] - tau[1,ja]
        dtau[2] = tau[2,ia] - tau[2,ja]
        dtau[3] = tau[3,ia] - tau[3,ja]
        
        isp = atm2species[ia]
        jsp = atm2species[ja]
        ZiZj = Zvals[isp]*Zvals[jsp]

        for i in -mmm1:mmm1, j in -mmm2:mmm2, k in -mmm3:mmm3
            T[1] = i*t1[1] + j*t2[1] + k*t3[1]
            T[2] = i*t1[2] + j*t2[2] + k*t3[2]
            T[3] = i*t1[3] + j*t2[3] + k*t3[3]
            D[1] = dtau[1] - T[1]
            D[2] = dtau[2] - T[2]
            D[3] = dtau[3] - T[3]
            D2 = D[1]^2 + D[2]^2 + D[3]^2
            Dmag = sqrt(D2)
            xx = _Herfc(η*Dmag)
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

    F_NN[:] = F_NN_G[:] + F_NN_R[:]

    return
end