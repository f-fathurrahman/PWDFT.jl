import PWDFT: Poisson_solve

#
# Given electron density in real space, return Hartree potential in reciprocal
# space
#
function Poisson_solve( pw::PWGridGamma, rhoR::Array{Float64,1} )
    gvec = pw.gvec
    G2 = gvec.G2
    Ng = gvec.Ng
    idx_g2r = gvec.idx_g2r
    idx_g2rm = gvec.idx_g2rm

    ctmp = R_to_G( pw, rhoR )
    
    ctmp[1] = 0.0 + im*0.0  # the first GVectors is zero vector

    for ig = 2:Ng
        #
        ip = idx_g2r[ig]
        ctmp[ip] = ctmp[ip]/G2[ig]
        #
        ipm = idx_g2rm[ig]
        ctmp[ipm] = conj(ctmp[ipm])/G2[ig]
    end
    return 4.0*pi*ctmp
end