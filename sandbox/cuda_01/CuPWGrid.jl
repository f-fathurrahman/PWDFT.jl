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

    return CuGVectorsW( Ngwx, Ngw, idx_gw2g, idx_gw2g, kpoints )

end

# CuArray version of PWGrid
struct CuPWGrid
    ecutwfc::Float64
    ecutrho::Float64
    Ns::Tuple{Int64,Int64,Int64}
    LatVecs::Array{Float64,2}
    RecVecs::Array{Float64,2}
    CellVolume::Float64
    r::Array{Float64,2} # not really used for now
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
    r = PWDFT.init_grid_R( Ns, LatVecs )  # FIXME: Not really used for now
    
    gvec_ = PWDFT.init_gvec( Ns, RecVecs, ecutrho )
    gvec = CuGVectors( gvec_ )

    if kpoints == nothing
        kpoints = KPoints( 1, (1,1,1), zeros(3,1), [1.0], RecVecs )
    end

    gvecw_ = PWDFT.init_gvecw( ecutwfc, gvec_, kpoints )
    gvecw = CuGVectorsW( gvecw_ )

    planfw = plan_fft( CuArrays.fill(0.0 + im*0.0, Ns) )
    planbw = plan_ifft( CuArrays.fill(0.0 + im*0.0, Ns) )

    return CuPWGrid( ecutwfc, ecutrho, Ns, LatVecs, RecVecs, CellVolume, r, gvec, gvecw,
                     planfw, planbw )
end
