using CuArrays

using PWDFT

struct CuGVectors
    Ng::Int64
    G::CuArray{Float64,2}
    G2::CuArray{Float64,1}
    idx_g2r::CuArray{Int64,1}
    G2_shells::CuArray{Float64,1}
    idx_g2shells::CuArray{Int64,1}
end

function CuGVectors( gvec::GVectors )

    G = CuArray( gvec.G )
    G2 = CuArray( gvec.G2 )
    idx_g2r = CuArray( gvec.idx_g2r )
    G2_shells = CuArray( gvec.G2_shells )
    idx_g2shells  = CuArray( gvec.idx_g2shells )

    return CuGVectors(gvec.Ng, G, G2, idx_g2r, G2_shells, idx_g2shells)
end




struct CuGVectorsW
    Ngwx::Int64
    Ngw::Array{Int64,1}
    idx_gw2g::Array{CuArray{Int64,1},1}
    idx_gw2r::Array{CuArray{Int64,1},1}
    kpoints::KPoints
end

"""
Initialize CuGVectorsW from an instance of `GVectorsW` by direct copy.
"""
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

function main()

    pw = PWGrid(15.0, gen_lattice_fcc(10.0))
    
    cu_gvec = CuGVectors( pw.gvec )
    cu_gvecw = CuGVectorsW( pw.gvecw )

    println("Pass here")
end

@time main()
@time main()