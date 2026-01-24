# This is assuming that rho1 and rho2 are in (RhoeTotal, magn)

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

    fac = 4π

    idx_g2r = pw.gvec.idx_g2r
    G2 = pw.gvec.G2

    Nspin = size(rho1, 2)
    @assert Nspin == size(rho2, 2)
    
    res = 0.0
    for ig in 2:Ngf
        ip = idx_g2r[ig]
        res += real( conj(rho1[ip,1]) * rho2[ip,1] ) / G2[ig]
    end
    res *= fac

    # XXX This also includes the case of noncollinear magn
    if Nspin >= 2
        res += fac*sum(real.( conj(rho1[1,2:Nspin]) .* rho2[1,2:Nspin] )) # why?
        for ig in 2:Ngf
            ip = idx_g2r[ig]
            res += fac*sum(real.( conj(rho1[ip,2:Nspin]) .* rho2[ip,2:Nspin] ))
        end
    end
    return 0.5 * res * pw.CellVolume # XXX need factor of 1/2 ?
end


function _sum_until_Ngf(pw, r::Array{Float64})
    if pw.using_dual_grid
        Ngf = pw.gvecs.Ng
    else
        Ngf = pw.gvec.Ng
    end
    idx_g2r = pw.gvec.idx_g2r
    ss = 0.0 + im*0.0
    for ig in 1:Ngf
        ip = idx_g2r[ig]
        ss += r[ip]
    end
    return ss
end