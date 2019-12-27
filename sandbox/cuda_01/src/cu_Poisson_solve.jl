import PWDFT: Poisson_solve

function kernel_Poisson_solve!( idx_g2r, G2, ctmp )
    
    ig = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x

    Ng = length(G2)

    if (ig <= Ng) && (ig != 1)
        ip = idx_g2r[ig]
        ctmp[ip] = ctmp[ip]/G2[ig]
    end
    
    return
end


function Poisson_solve( pw::CuPWGrid, RhoeR::CuArray{Float64,1} )

    G2 = pw.gvec.G2
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r
    Npoints = prod(pw.Ns)
    
    ctmp = CuArrays.zeros( ComplexF64, Npoints )
    ctmp[2:Npoints] = R_to_G( pw, RhoeR )[2:Npoints]
    
    #ctmp = R_to_G( pw, RhoeR )
    #@allowscalar ctmp[1] = 0.0 + 0.0*im

    Nthreads = 256
    Nblocks = ceil( Int64, Ng/Nthreads )

    @cuda threads=Nthreads blocks=Nblocks kernel_Poisson_solve!( idx_g2r, G2, ctmp )

    return 4.0*pi*ctmp # this is phiG
end