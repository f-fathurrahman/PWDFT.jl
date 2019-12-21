import PWDFT: op_K

# Factor of 0.5 is included
function kernel_op_K!( ist, G, idx_gw2g_k, k1, k2, k3, psi, Kpsi )

    igk = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x

    if igk <= length(idx_gw2g_k)
        ig = idx_gw2g_k[igk]
        Gw2 = (G[1,ig] + k1)^2 + (G[2,ig] + k2)^2 + (G[3,ig] + k3)^2
        Kpsi[igk,ist] = 0.5*psi[igk,ist]*Gw2
    end

    return
end

function op_K( Ham::CuHamiltonian, psi::CuArray{ComplexF64,2} )

    # Get global index of current k-point index
    ik = Ham.ik

    Nstates = size(psi)[2]

    Ngw = Ham.pw.gvecw.Ngw
    idx_gw2g_k = Ham.pw.gvecw.idx_gw2g[ik]
    G = Ham.pw.gvec.G
    k1 = Ham.pw.gvecw.kpoints.k[1,ik]
    k2 = Ham.pw.gvecw.kpoints.k[2,ik]
    k3 = Ham.pw.gvecw.kpoints.k[3,ik]

    Kpsi = CuArrays.zeros( ComplexF64, size(psi) )

    Nthreads = 256
    Nblocks = ceil( Int64, Ngw[ik]/Nthreads )

    for ist = 1:Nstates
        @cuda threads=Nthreads blocks=Nblocks kernel_op_K!( ist, G, idx_gw2g_k, k1, k2, k3, psi, Kpsi )
    end

    return Kpsi # factor of 0.5 is already included
end


function op_K( Ham::CuHamiltonian, psiks::CuBlochWavefunc )
    Nstates = size(psiks[1])[2] # Nstates should be similar for all Bloch states
    
    Nspin = Ham.electrons.Nspin
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    out = zeros_CuBlochWavefunc( Ham )
    
    for ispin = 1:Nspin, ik = 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        out[ikspin] = op_K( Ham, psiks[ikspin] )
    end
    return out
end
