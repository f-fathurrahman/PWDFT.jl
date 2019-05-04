function gen_V_Ps_loc_short( Ham::Hamiltonian )

    Nspecies = Ham.atoms.Nspecies
    Npoints = prod(Ham.pw.Ns)
    CellVolume = Ham.pw.CellVolume
    G2 = Ham.pw.gvec.G2

    V_Ps_loc_short = zeros(Npoints)
    Vg = zeros(ComplexF64, Npoints)

    strf = calc_strfact( Ham.atoms, Ham.pw )

    for isp = 1:Nspecies
        pspot = Ham.pspots[isp]
        for ig = 1:Ham.pw.gvec.Ng
            ip = Ham.pw.gvec.idx_g2r[ig]
            Vg[ip] = strf[ig,isp] * PWDFT.eval_Vloc_G_short( pspot, G2[ig] )
        end
        #
        V_Ps_loc_short[:] = V_Ps_loc_short[:] + real( G_to_R(Ham.pw, Vg) ) * Npoints / CellVolume
    end

    return V_Ps_loc_short

end


function gen_Rhoe_aux_G( Ham:: Hamiltonian )
    return gen_Rhoe_aux_G( Ham.atoms, Ham.pw )
end

function gen_Rhoe_aux_G( atoms::Atoms, pw::PWGrid, pspots; TOL = 1e-8 )

    Zvals = atoms.Zvals

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

    for isp = 1:Nspecies
        rlocal = pspots[isp].rlocal
        alpha = 1/(rlocal*sqrt(2))
        eta = sqrt(0.5*alpha)
        beta_s = rlocal*sqrt(2)
        for ig = 1:Ng
            ip = idx_g2r[ig]
            # determine eta
            #eta = 0.5*Gcut^2/-log(TOL)
            #alpha = 2*eta^2  # alternative parameter
            #rloc = 1/(alpha*sqrt(2))
            #Rhoe_aux_G[ip] = Rhoe_aux_G[ip] - Zvals[isp]*exp(-0.125*G2[ig]/eta^2)*Sf[ig,isp]/CellVolume
            Rhoe_aux_G[ip] = Rhoe_aux_G[ip] - Zvals[isp]*exp( -0.25*G2[ig]*beta_s^2 )*Sf[ig,isp]/CellVolume
        end
    end

    Rhoe_aux = real( G_to_R(pw,Rhoe_aux_G) )

    println("Rhoe_aux_G(G=0) = ", Rhoe_aux_G[1]*CellVolume)

    E_self = 0.0
    for ia = 1:atoms.Natoms
        isp = atoms.atm2species[ia]
        rlocal = pspots[isp].rlocal
        #alpha = 1/(rlocal*sqrt(2))
        beta_s = rlocal*sqrt(2)
        #E_self = E_self + sqrt(alpha/pi)*Zvals[isp]^2
        E_self = E_self + 1/sqrt(2*pi)/beta_s*Zvals[isp]^2  # beta_s = 1
    end
    
    return Rhoe_aux*Npoints, E_self #/CellVolume
end



function gen_V_aux_G( Ham )
    return gen_V_aux_G( Ham.atoms, Ham.pw, Ham.pspots )
end

function gen_V_aux_G( atoms::Atoms, pw::PWGrid, pspots; TOL = 1e-8 )

    Zvals = atoms.Zvals

    # determine eta
    Gcut = 2*pw.ecutwfc/(2*pi)
    eta = 0.5*Gcut^2/-log(TOL)

    alpha = 2*eta^2

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

    for isp = 1:Nspecies
        
        psp = pspots[isp]
        rlocal = psp.rlocal
        beta_s = rlocal*sqrt(2)

        for ig = 1:Ng

            ip = idx_g2r[ig]

            Gm = sqrt(G2[ig])
            Gr = Gm*rlocal

            expGr2 = exp(-0.5*Gr^2)      
            if Gm > eps()
                V_aux_G[ip] = V_aux_G[ip] - 4*pi*Zvals[isp]*expGr2/G2[ig]*Sf[ig,isp]/CellVolume
            else
                #myfunc(r) = r^2 * ( PWDFT.eval_Vloc_R(psp, r) + Zvals[isp]*erf(r/beta_s) )
                #V0 = quadgk( myfunc, eps(), 10.0 )[1]
                #V_aux_G[ip] = V_aux_G[ip] - 4*pi*V0*Sf[ig,isp]/CellVolume
                #V_aux_G[ip] = V_aux_G[ip] - 2*pi*Zvals[isp]*rloc^2*Sf[ig,isp]/CellVolume
                #V_aux_G[ip] = V_aux_G[ip] - 2*pi*Zvals[isp]*rloc^2*Sf[ig,isp]/CellVolume
                #V_aux_G[ip] = V_aux_G[ip] - 4*pi*(2*Zvals[isp]/sqrt(pi))*Sf[ig,isp]/CellVolume
            end
        end
    end
    println("V_aux_G[1] = ", V_aux_G[1])
    println("V_aux_G[2] = ", V_aux_G[2])

    V_aux = real(G_to_R(pw, V_aux_G))*Npoints

    println("integ V_aux in R-space, constructed from G-space: ", sum(V_aux)*dVol)

    return V_aux  # note the minus sign
end



function calc_E_alphat( atoms::Atoms, pw::PWGrid, pspots; TOL = 1e-8 )

    Zvals = atoms.Zvals
    atm2species = atoms.atm2species
    Natoms = atoms.Natoms

    # determine eta
    Gcut = 2*pw.ecutwfc/(2*pi)
    eta = 0.5*Gcut^2/-log(TOL)
    alpha = 2*eta^2

    # structure factor
    Sf = calc_strfact( atoms, pw )
    
    Npoints = prod(pw.Ns)
    CellVolume = pw.CellVolume
    dVol = CellVolume/Npoints

    E_alphat = 0.0
    for ia = 1:Natoms
        isp = atm2species[ia]
        psp = pspots[isp]
        rlocal = psp.rlocal
        beta_s = rlocal*sqrt(2)

        myfunc(r) = r^2 * ( PWDFT.eval_Vloc_R(psp, r) + Zvals[isp]*erf(r/beta_s) )
        #myfunc(r) = r^2 * ( PWDFT.eval_Vloc_R(psp, r) + Zvals[isp]/r )
        
        E_alphat = E_alphat + 4*pi*quadgk( myfunc, eps(), 10.0 )[1]*Zvals[isp]
    end

    return E_alphat/CellVolume/Npoints
end




function gen_V_Ps_loc_screened( atoms::Atoms, pw::PWGrid, pspots; TOL = 1e-8 )

    Zvals = atoms.Zvals

    # determine eta
    Gcut = 2*pw.ecutwfc/(2*pi)
    eta = 0.5*Gcut^2/-log(TOL)

    alpha = 2*eta^2

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

    println("Generating V_Ps_loc screened")

    for isp = 1:Nspecies
        
        psp = pspots[isp]
        rlocal = psp.rlocal
        beta_s = rlocal*sqrt(2)

        for ig = 1:Ng

            ip = idx_g2r[ig]
            
            Gm = sqrt(G2[ig])

            myfunc(r) = r^2 * besselj0(r*Gm) * ( PWDFT.eval_Vloc_R(psp, r) + Zvals[isp]*erf(r/beta_s) )
            
            V0 = quadgk( myfunc, eps(), 20.0 )[1]

            #@printf("Gm = %18.10f, V0 = %18.10f\n", Gm, V0)

            V_aux_G[ip] = V_aux_G[ip] + 4*pi*V0*Sf[ig,isp]/CellVolume
        end
    end

    V_aux = real(G_to_R(pw, V_aux_G))*Npoints

    return V_aux
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

function gen_Rhoe_aux_R( Ham::Hamiltonian )
    return gen_Rhoe_aux_R( Ham.atoms, Ham.pw )
end

function gen_V_aux_R( Ham::Hamiltonian )
    return gen_V_aux_R( Ham.atoms, Ham.pw )
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


function gen_Rhoe_aux_G( atoms::Atoms, pw::PWGrid; TOL = 1e-8 )

    Zvals = atoms.Zvals

    # determine eta
    Gcut = 2*pw.ecutwfc/(2*pi)
    eta = 0.5*Gcut^2/-log(TOL)

    alpha = 2*eta^2  # alternative parameter

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

    ss = 0.0
    for ia = 1:atoms.Natoms
        isp = atoms.atm2species[ia]
        ss = ss + Zvals[isp]^2
    end
    E_self = sqrt(alpha/pi)*ss
    
    return Rhoe_aux*Npoints, E_self/CellVolume
end
