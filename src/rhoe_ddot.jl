function rhoe_ddot( pw, rho1, rho2 )
#=
  !! Calculates \(4\pi/G^2\ \rho_1(-G)\ \rho_2(G) = V1_\text{Hartree}(-G)\ \rho_2(G)\)
  !! used as an estimate of the self-consistency error on the energy.
=#

    if pw.using_dual_grid
        Ngf = pw.gvecs.Ng
    else
        Ngf = pw.gvec.Ng
    end

    @info "rhoe_ddot: Ngf = $(Ngf)"

    fac = 4π

    idx_g2r = pw.gvec.idx_g2r
    G2 = pw.gvec.G2

    Nspin = size(rho1, 2)
    @assert Nspin == size(rho2, 2)

    # Total RhoeG

    if Nspin == 2
        ρ1 = rho1[:,1] + rho1[:,2]
        ρ2 = rho2[:,1] + rho2[:,2]
    else
        ρ1 = rho1[:,1]
        ρ2 = rho2[:,1]
    end

    @info "rhoe_ddot: sum ρ1 = $(_sum_until_Ngf(pw, ρ1))"
    @info "rhoe_ddot: sum ρ2 = $(_sum_until_Ngf(pw, ρ2))"
    
    res = 0.0
    for ig in 2:Ngf
        ip = idx_g2r[ig]
        res += real( conj(ρ1[ip,1]) * ρ2[ip,1] ) / G2[ig]
    end
    res *= fac
    println("res in rhoe_ddot after total rhoe = ", res)

    # This is the magnetization?
    if Nspin == 2
        μ1 = rho1[:,1] - rho1[:,2]
        μ2 = rho2[:,1] - rho2[:,2]
        @info "rhoe_ddot: sum μ1 = $(_sum_until_Ngf(pw, μ1))"
        @info "rhoe_ddot: sum μ2 = $(_sum_until_Ngf(pw, μ2))"
        #ip = 1 # ig = 1
        res += fac*real( conj(μ1[1]) * μ2[1] ) # why?
        @info "res line 52 = $(res)"
        for ig in 2:Ngf
            ip = idx_g2r[ig]
            res += fac*real( conj(μ1[ip]) * μ2[ip] )
        end
        @info "res line 57 = $(res)"
    end
    println("res in rhoe_ddot after magn = ", 0.5 * res * pw.CellVolume)
    return 0.5 * res * pw.CellVolume # XXX need factor of 1/2 ?
end


function _sum_until_Ngf(pw, ρ)
    if pw.using_dual_grid
        Ngf = pw.gvecs.Ng
    else
        Ngf = pw.gvec.Ng
    end
    idx_g2r = pw.gvec.idx_g2r
    ss = 0.0 + im*0.0
    for ig in 1:Ngf
        ip = idx_g2r[ig]
        ss += ρ[ip]
    end
    return ss
end