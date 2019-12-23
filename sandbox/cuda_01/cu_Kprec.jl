import PWDFT: Kprec


function kernel_Kprec!( ist, idx_gw2g_ik, G, k1, k2, k3, psi, Kpsi )
    
    igk = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x

    Ngw_ik = length( idx_gw2g_ik )

    if igk <= Ngw_ik
        
        ig = idx_gw2g_ik[igk]
        
        Gw2 = ( G[1,ig] + k1 )^2 + ( G[2,ig] + k2 )^2 + ( G[3,ig] + k3 )^2
            
        Kpsi[igk,ist] = psi[igk,ist] / ( 1.0 + Gw2 )
    
    end
    
    return
end


function Kprec( ik::Int64, pw::CuPWGrid, psi )

    Ngw_ik  = size(psi)[1]
    Nstates = size(psi)[2]
    idx_gw2g_ik = pw.gvecw.idx_gw2g[ik]
    G = pw.gvec.G
    k1 = pw.gvecw.kpoints.k[1,ik]
    k2 = pw.gvecw.kpoints.k[2,ik]
    k3 = pw.gvecw.kpoints.k[3,ik]

    Kpsi = CuArrays.zeros( ComplexF64, Ngw_ik, Nstates )

    Nthreads = 256
    Nblocks = ceil(Int64, Ngw_ik/Nthreads)
    for ist = 1:Nstates
        @cuda threads=Nthreads blocks=Nblocks kernel_Kprec!( ist, idx_gw2g_ik, G, k1, k2, k3, psi, Kpsi )
    end
    return Kpsi
end

