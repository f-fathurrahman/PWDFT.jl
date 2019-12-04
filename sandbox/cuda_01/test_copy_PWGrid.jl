using CuArrays

using PWDFT

struct CuGVectors
    Ng::Int64
    G::Array{Float64,2}
    G2::Array{Float64,1}
    idx_g2r::Array{Int64,1}
    G2_shells::Array{Float64,1}
    idx_g2shells::Array{Int64,1}
end

struct CuGVectorsW
    Ngwx::Int64
    Ngw::Array{Int64,1}
    idx_gw2g  # two dim array the size if idx_gw2r
    idx_gw2r
    kpoints
end


function init_CuGVectorsW( gvecw::GVectorsW )
    
    Ngwx = gvecw.Ngwx
    Ngw = gvecw.Ngw

    kpoints = gvecw.kpoints
    Nkpt = kpoints.Nkpt
    
    idx_gw2g = CuArrays.fill(0, (Ngwx, Nkpt))
    idx_gw2r = CuArrays.fill(0, (Ngwx, Nkpt))
    
    # Copy to GPU
    for ik in 1:Nkpt
        Ngwk = Ngw[ik]
        idx_gw2g[1:Ngwk] = CuArray( gvecw.idx_gw2g[ik] )
        idx_gw2r[1:Ngwk] = CuArray( gvecw.idx_gw2r[ik] )
    end
    
    println(typeof(idx_gw2g))
    println(typeof(idx_gw2r))

    return CuGVectorsW( Ngwx, Ngw, idx_gw2g, idx_gw2g, kpoints )

end

function main()

    pw = PWGrid(15.0, gen_lattice_fcc(10.0))
    init_CuGVectorsW( pw.gvecw)

    println("Pass here")
end

@time main()
@time main()