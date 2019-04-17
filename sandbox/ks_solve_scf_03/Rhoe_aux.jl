function gen_Rhoe_aux_G( Ham:: Hamiltonian )
    return gen_Rhoe_aux_G( Ham.atoms, Ham.pw )
end

function gen_Rhoe_aux_G( atoms::Atoms, pw::PWGrid; TOL = 1e-8 )

    Zvals = atoms.Zvals

    # determine eta
    Gcut = 2*pw.ecutwfc/(2*pi)
    eta = 0.5*Gcut^2/-log(TOL)

    # structure factor
    Sf = calc_strfact( atoms, pw )
    
    Ng = pw.gvec.Ng
    G2 = pw.gvec.G2
    Nspecies = atoms.Nspecies
    Npoints = prod(pw.Ns)
    idx_g2r = pw.gvec.idx_g2r
    CellVolume = pw.CellVolume
    dVol = CellVolume/Npoints

    Rhoe_aux_G = zeros(ComplexF64,Npoints)
    Rhoe_aux = zeros(Float64,Npoints)

    for ig = 1:Ng
        ip = idx_g2r[ig]
        for isp = 1:Nspecies
            Rhoe_aux_G[ip] = Rhoe_aux_G[ip] - Zvals[isp]*exp(-0.125*G2[ig]/eta^2)*Sf[ig,isp]/CellVolume
        end
    end

    Rhoe_aux = real( G_to_R(pw,Rhoe_aux_G) )

    println("Rhoe_aux_G(G=0) = ", Rhoe_aux_G[1]*CellVolume)
    
    return Rhoe_aux*Npoints  # note the minus sign
end


function gen_Rhoe_aux_R( Ham::Hamiltonian )
    return gen_Rhoe_aux_R( Ham.atoms, Ham.pw )
end

function gen_V_aux_R( Ham::Hamiltonian )
    return gen_V_aux_R( Ham.atoms, Ham.pw )
end



# only intended for testing
function gen_Rhoe_aux_R( atoms::Atoms, pw::PWGrid; TOL = 1e-8 )

    Zvals = atoms.Zvals

    # determine eta
    Gcut = 2*pw.ecutwfc/(2*pi)
    eta = 0.5*Gcut^2/-log(TOL)

    alpha = 2*eta^2  # alternative parameter
    
    Nspecies = atoms.Nspecies
    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume
    dVol = CellVolume/Npoints
    Natoms = atoms.Natoms

    Rhoe_aux = zeros(Float64,Npoints)
    for ip = 1:Npoints
        r = pw.r[:,ip]
        for ia = 1:Natoms
            isp = atoms.atm2species[ia]
            R = atoms.positions[:,ia]
            dr2 = dot(r-R,r-R)
            #Rhoe_aux[ip] = Rhoe_aux[ip] + Zvals[isp]*exp(-2*eta^2*dr2)
            Rhoe_aux[ip] = Rhoe_aux[ip] + Zvals[isp]*exp(-alpha*dr2)
        end
    end
    
    #Rhoe_aux = -(2*eta)^3/((2*pi)^1.5)*Rhoe_aux
    Rhoe_aux = -(alpha/pi)^(3/2) * Rhoe_aux
    
    println("integ Rhoe_aux in R-space: ", sum(Rhoe_aux)*dVol)
    return Rhoe_aux
end


# only intended for testing
function gen_V_aux_R( atoms::Atoms, pw::PWGrid; TOL = 1e-8 )

    Zvals = atoms.Zvals

    # determine eta
    Gcut = 2*pw.ecutwfc/(2*pi)
    eta = 0.5*Gcut^2/-log(TOL)

    alpha = 2*eta^2

    Nspecies = atoms.Nspecies
    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume
    dVol = CellVolume/Npoints
    Natoms = atoms.Natoms

    t1m = dot(atoms.LatVecs[:,1],atoms.LatVecs[:,1])
    t2m = dot(atoms.LatVecs[:,2],atoms.LatVecs[:,2])
    t3m = dot(atoms.LatVecs[:,3],atoms.LatVecs[:,3])

    tmax = sqrt(0.5*(-log(TOL)))/eta
    mmm1 = round(Int64, tmax/t1m + 1.5) + 1
    mmm2 = round(Int64, tmax/t2m + 1.5) + 1
    mmm3 = round(Int64, tmax/t3m + 1.5) + 1

    V_aux = zeros(Float64,Npoints)

    t1 = atoms.LatVecs[:,1]
    t2 = atoms.LatVecs[:,2]
    t3 = atoms.LatVecs[:,3]
    T = zeros(3)
    R = zeros(3)

    @printf("Real space summation: %d %d %d\n", mmm1, mmm2, mmm3)

    for i = -mmm1:mmm1
    for j = -mmm2:mmm2
    for k = -mmm3:mmm3
        T[1] = i*t1[1] + j*t2[1] + k*t3[1]
        T[2] = i*t1[2] + j*t2[2] + k*t3[2]
        T[3] = i*t1[3] + j*t2[3] + k*t3[3]
        for ip = 1:Npoints
            r = pw.r[:,ip]
            for ia = 1:Natoms
                isp = atoms.atm2species[ia]
                R = atoms.positions[:,ia] - T
                dr = sqrt(dot(r-R,r-R))
                if dr > eps()
                    V_aux[ip] = V_aux[ip] + Zvals[isp]*erf(sqrt(alpha)*dr)/dr
                else
                    println("Pass here")
                    V_aux[ip] = V_aux[ip] + 2*Zvals[isp]*sqrt(alpha/pi)
                end
            end
        end
    end
    end
    end
    
    println("integ V_aux in R-space: ", sum(V_aux)*dVol)
    return V_aux
end




function gen_V_aux_G( atoms::Atoms, pw::PWGrid; TOL = 1e-8 )

    Zvals = atoms.Zvals

    # determine eta
    Gcut = 2*pw.ecutwfc/(2*pi)
    eta = 0.5*Gcut^2/-log(TOL)

    alpha = 2*eta^2

    #alpha = 1/(rloc*sqrt(2))
    rloc = 1/(alpha*sqrt(2))

    # structure factor
    Sf = calc_strfact( atoms, pw )
    
    Ng = pw.gvec.Ng
    G2 = pw.gvec.G2
    Nspecies = atoms.Nspecies
    Npoints = prod(pw.Ns)
    idx_g2r = pw.gvec.idx_g2r
    CellVolume = pw.CellVolume
    dVol = CellVolume/Npoints

    V_aux_G = zeros(ComplexF64,Npoints)
    V_aux = zeros(Float64,Npoints)


    for ig = 2:Ng
        ip = idx_g2r[ig]

        Gm = sqrt(G2[ig])
        Gr = Gm*rloc
        expGr2 = exp(-0.5*Gr^2)

        for isp = 1:Nspecies
            if Gm > eps()
                V_aux_G[ip] = V_aux_G[ip] - 4*pi*Zvals[isp]*expGr2/G2[ig]*Sf[ig,isp]/CellVolume
            else
                V_aux_G[ip] = V_aux_G[ip] + 2*pi*Zvals[isp]*rloc^2*Sf[ig,isp]/CellVolume
            end
        end
    end


    V_aux = G_to_R(pw, V_aux_G)*Npoints

    println("integ V_aux in R-space, constructed from G-space: ", sum(V_aux)*dVol)

    return V_aux  # note the minus sign
end
