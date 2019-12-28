# CuArray version of GVectors
struct CuGVectors
    Ng::Int64
    G::CuArray{Float64,2}
    G2::CuArray{Float64,1}
    idx_g2r::CuArray{Int64,1}
    G2_shells::CuArray{Float64,1}
    idx_g2shells::CuArray{Int64,1}
end

#
# Simple copy version
#
function CuGVectors( gvec::GVectors )
    G = CuArray( gvec.G )
    G2 = CuArray( gvec.G2 )
    idx_g2r = CuArray( gvec.idx_g2r )
    G2_shells = CuArray( gvec.G2_shells )
    idx_g2shells = CuArray( gvec.idx_g2shells )

    return CuGVectors(gvec.Ng, G, G2, idx_g2r, G2_shells, idx_g2shells)
end


# CuArray version of GVectorsW
struct CuGVectorsW
    Ngwx::Int64
    Ngw::Array{Int64,1}
    idx_gw2g::Array{CuArray{Int64,1},1}
    idx_gw2r::Array{CuArray{Int64,1},1}
    kpoints::KPoints
end


# Naive copy
function CuGVectorsW( gvecw::GVectorsW )

    Ngwx = gvecw.Ngwx
    Ngw = gvecw.Ngw

    kpoints = gvecw.kpoints
    Nkpt = kpoints.Nkpt

    TYP = typeof(CuArray([1]))
    idx_gw2g = Array{TYP,1}(undef,Nkpt)
    idx_gw2r = Array{TYP,1}(undef,Nkpt)

    # Copy to GPU
    for ik in 1:Nkpt
        idx_gw2g[ik] = CuArray( gvecw.idx_gw2g[ik] )
        idx_gw2r[ik] = CuArray( gvecw.idx_gw2r[ik] )
    end

    return CuGVectorsW( Ngwx, Ngw, idx_gw2g, idx_gw2r, kpoints )

end

# CuArray version of PWGrid
struct CuPWGrid
    ecutwfc::Float64
    ecutrho::Float64
    Ns::Tuple{Int64,Int64,Int64}
    LatVecs::Array{Float64,2}
    RecVecs::Array{Float64,2}
    CellVolume::Float64
    gvec::CuGVectors
    gvecw::CuGVectorsW
    planfw::CuArrays.CUFFT.cCuFFTPlan{Complex{Float64},-1,false,3}
    planbw::AbstractFFTs.ScaledPlan{Complex{Float64},CuArrays.CUFFT.cCuFFTPlan{Complex{Float64},1,false,3},Float64}
end


function CuPWGrid( ecutwfc::Float64, LatVecs::Array{Float64,2}; kpoints=nothing, Ns_=(0,0,0) )

    ecutrho = 4.0*ecutwfc

    RecVecs = 2*pi*inv(Matrix(LatVecs'))

    CellVolume = abs(det(LatVecs))
    #
    LatVecsLen = Array{Float64}(undef,3)
    LatVecsLen[1] = norm(LatVecs[:,1])
    LatVecsLen[2] = norm(LatVecs[:,2])
    LatVecsLen[3] = norm(LatVecs[:,3])

    Ns1 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[1]/pi ) + 1
    Ns2 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[2]/pi ) + 1
    Ns3 = 2*round( Int64, sqrt(ecutrho/2)*LatVecsLen[3]/pi ) + 1

    if any(Ns_ .== 0)
        Ns1 = good_fft_order(Ns1)
        Ns2 = good_fft_order(Ns2)
        Ns3 = good_fft_order(Ns3)
        Ns = (Ns1,Ns2,Ns3)
    else
        Ns = Ns_[:]
    end

    Npoints = prod(Ns)

    gvec_ = PWDFT.init_gvec( Ns, RecVecs, ecutrho )
    gvec = CuGVectors( gvec_ )

    if kpoints == nothing
        kpoints = KPoints( 1, (1,1,1), zeros(3,1), [1.0], RecVecs )
    end

    gvecw_ = PWDFT.init_gvecw( ecutwfc, gvec_, kpoints )
    gvecw = CuGVectorsW( gvecw_ )

    planfw = plan_fft( CuArrays.zeros(ComplexF64, Ns) )
    planbw = plan_ifft( CuArrays.zeros(ComplexF64, Ns) )

    return CuPWGrid( ecutwfc, ecutrho, Ns, LatVecs, RecVecs, CellVolume, gvec, gvecw,
                     planfw, planbw )
end

import PWDFT: op_nabla, op_nabla_dot

function kernel_op_nabla!( idx_g2r, G, RhoeG, res )
    ig = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x
    Ng = length(idx_g2r)
    if ig <= Ng
        ip = idx_g2r[ig]
        res[1,ip] = im*G[1,ig]*RhoeG[ip]
        res[2,ip] = im*G[2,ig]*RhoeG[ip]
        res[3,ip] = im*G[3,ig]*RhoeG[ip]
    end
    return
end

function op_nabla( pw::CuPWGrid, Rhoe::CuArray{Float64,1} )
    G = pw.gvec.G
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r
    Npoints = prod(pw.Ns)

    RhoeG = R_to_G(pw,Rhoe)  # We are not taking indexing here

    ∇RhoeG_full = CuArrays.zeros(ComplexF64,3,Npoints)
    ∇Rhoe = CuArrays.zeros(Float64,3,Npoints)

    Nthreads = 256
    Nblocks = ceil(Int64, Ng/Nthreads)

    @cuda threads=Nthreads blocks=Nblocks kernel_op_nabla!( idx_g2r, G, RhoeG, ∇RhoeG_full )

    ∇Rhoe[1,:] = real(G_to_R(pw,∇RhoeG_full[1,:]))
    ∇Rhoe[2,:] = real(G_to_R(pw,∇RhoeG_full[2,:]))
    ∇Rhoe[3,:] = real(G_to_R(pw,∇RhoeG_full[3,:]))
    return ∇Rhoe

end

function kernel_op_nabla_dot!( idx_g2r, G, hG, res )
    ig = ( blockIdx().x - 1 )*blockDim().x + threadIdx().x
    Ng = length(idx_g2r)
    if ig <= Ng
        ip = idx_g2r[ig]
        res[ip] = im*( G[1,ig]*hG[1,ip] + G[2,ig]*hG[2,ip] + G[3,ig]*hG[3,ip] )
    end
    return
end

function op_nabla_dot( pw::CuPWGrid, h::CuArray{Float64,2} )
    G = pw.gvec.G
    Ng = pw.gvec.Ng
    idx_g2r = pw.gvec.idx_g2r
    Npoints = prod(pw.Ns)

    # hG is using Npoints instead of Ng
    hG = CuArrays.zeros(ComplexF64,3,Npoints)
    hG[1,:] = R_to_G( pw, h[1,:] )
    hG[2,:] = R_to_G( pw, h[2,:] )
    hG[3,:] = R_to_G( pw, h[3,:] )

    divhG_full = CuArrays.zeros(ComplexF64,Npoints)

    Nthreads = 256
    Nblocks = ceil(Int64, Ng/Nthreads)

    @cuda threads=Nthreads blocks=Nblocks kernel_op_nabla_dot!( idx_g2r, G, hG, divhG_full )

    divh = real( G_to_R( pw, divhG_full ) )
    return divh

end
