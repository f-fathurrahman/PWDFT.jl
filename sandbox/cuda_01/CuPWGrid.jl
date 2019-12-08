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
    planfw::FFTW.cFFTWPlan{Complex{Float64},-1,false,3}
    planbw::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.cFFTWPlan{Complex{Float64},1,false,3},Float64}
end

