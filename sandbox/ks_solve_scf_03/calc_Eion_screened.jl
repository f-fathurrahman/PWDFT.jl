using SpecialFunctions: erfc

function calc_Eion_screened( atoms::Atoms )
    return calc_Eion_screened( atoms.LatVecs, atoms, atoms.Zvals )
end

function calc_Eion_screened( atoms::Atoms, Zvals::Array{Float64,1} )
    return calc_Eion_screened( atoms.LatVecs, atoms, Zvals )
end

function calc_Eion_screened( pw::PWGrid, atoms::Atoms, Zvals::Array{Float64,1} )
    return calc_Eion_screened( pw.LatVecs, atoms, Zvals )
end

function calc_Eion_screened( LatVecs::Array{Float64,2}, atoms::Atoms, Zvals::Array{Float64,1} )

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

    # Calculate η
    gcut = 2.0
    ebsl = 1e-8
    glast2 = gcut*gcut
    gexp = -log(ebsl)    
    η = sqrt(glast2/gexp)/2

    x = 0.0
    totalcharge = 0.0
    for ia = 1:Natoms
        isp = atm2species[ia]
        totalcharge = totalcharge + Zvals[isp]
    end

    tmax = sqrt(0.5*gexp)/η

    mmm1 = round(Int64, tmax/t1m + 1.5)
    mmm2 = round(Int64, tmax/t2m + 1.5)
    mmm3 = round(Int64, tmax/t3m + 1.5)

    dtau = zeros(Float64,3)
    T = zeros(Float64,3)

    Eion_scr = 0.0

    for ia = 1:Natoms
    for ja = 1:Natoms
        dtau[1] = tau[1,ia] - tau[1,ja]
        dtau[2] = tau[2,ia] - tau[2,ja]
        dtau[3] = tau[3,ia] - tau[3,ja]
        isp = atm2species[ia]
        jsp = atm2species[ja]
        ZiZj = Zvals[isp]*Zvals[jsp]
        for i = -mmm1:mmm1
        for j = -mmm2:mmm2
        for k = -mmm3:mmm3
            if (ia != ja) || ( (abs(i) + abs(j) + abs(k)) != 0 )
                T[1] = i*t1[1] + j*t2[1] + k*t3[1]
                T[2] = i*t1[2] + j*t2[2] + k*t3[2]
                T[3] = i*t1[3] + j*t2[3] + k*t3[3]
                rmag2 = sqrt( (dtau[1] - T[1])^2 +
                              (dtau[2] - T[2])^2 +
                              (dtau[3] - T[3])^2 )
                Eion_scr = Eion_scr + 0.5*ZiZj*erfc(rmag2*η)/rmag2
            end
        end
        end
        end
    end
    end

    return Eion_scr

end
