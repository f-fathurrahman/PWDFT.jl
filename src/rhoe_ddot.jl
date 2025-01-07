function rhoe_ddot( pw, rho1, rho2, Ngf )
#=
  !! Calculates \(4\pi/G^2\ \rho_1(-G)\ \rho_2(G) = V1_\text{Hartree}(-G)\ \rho_2(G)\)
  !! used as an estimate of the self-consistency error on the energy.
=#

    fac = 4π

    idx_g2r = pw.gvec.idx_g2r
    G2 = pw.gvec.G2

    # Total RhoeG
    res = 0.0
    if Nspin == 2
        ρ1 = rho1[:,1] + rho1[:,2]
        ρ2 = rho2[:,1] + rho2[:,2]
    else
        ρ1 = rho[:,1]
        ρ2 = rho[:,1]
    end
    for ig in 2:Ngf
        ip = idx_g2r[ig]
        res += real( conj(ρ1[ip,1]) * ρ2[ig,1] ) / G2[ig]
    end
    res *= fac
    
    # This is the magnetization?
    if Nspin == 2
        μ1 = rho1[:,1] - rho1[:,2]
        μ2 = rho2[:,1] - rho2[:,2]        
        #ip = 1 # ig = 1
        res += real( conj(μ1[1]) * μ2[1] ) # why?
        for ig in 2:Ngf
            ip = idx_g2r[ig]
            res += real( conj(μ1[ip]) * μ2[ip] )
        end
        res *= fac
    end
  
    return 0.5 * res * CellVolume # XXX need factor of 1/2 ?
end
