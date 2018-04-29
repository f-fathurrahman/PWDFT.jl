"""
Bloch wave vector
"""
mutable struct KPoints
    Nkpt::Int64
    k::Array{Float64,2}
    wk::Array{Float64,1}
    RecVecs::Array{Float64,2} # copy of reciprocal vectors
end


import Base.println

function println( kpoints::KPoints )

    @printf("\n")
    @printf("Total number of kpoints = %d\n", kpoints.Nkpt )
    ss = abs(kpoints.RecVecs[1,1])
    
    @printf("\n")
    @printf("kpoints in Cartesian coordinate (scale: %f)\n", ss)
    @printf("\n")
    kcart = kpoints.RecVecs*kpoints.k / ss
    for ik = 1:kpoints.Nkpt
        @printf("%4d [%14.10f %14.10f %14.10f] %14.10f\n",
                ik, kcart[1,ik], kcart[2,ik], kcart[3,ik], kpoints.wk[ik])
    end

    @printf("\n")
    @printf("kpoints in Cartesian coordinate (unscaled)\n")
    @printf("\n")
    kcart = kpoints.RecVecs*kpoints.k
    for ik = 1:kpoints.Nkpt
        @printf("%4d [%14.10f %14.10f %14.10f] %14.10f\n",
                ik, kcart[1,ik], kcart[2,ik], kcart[3,ik], kpoints.wk[ik])
    end
    @printf("\n")
    @printf("sum wk = %f\n", sum(kpoints.wk))
end


"""
Generate uniform kpoints-grid
"""
function KPoints( atoms::Atoms, mesh::Array{Int64,1}, is_shift::Array{Int64,1};
                  time_reversal=1, verbose=false )
    
    if verbose
        @printf("\n")
        @printf("Generating kpoints:\n")
        @printf("mesh     = (%d,%d,%d)\n", mesh[1], mesh[2], mesh[3])
        @printf("is_shift = (%d,%d,%d)\n", is_shift[1], is_shift[2], is_shift[3])
    end

    Nkpt, kgrid, mapping =
    spg_get_ir_reciprocal_mesh( atoms, mesh, is_shift, is_time_reversal=time_reversal )
    
    # search for unique mapping
    umap = unique(mapping)

    # Total number of grid points (unreduced)
    NkptTotal = prod(mesh)

    list_ir_k = []
    for ikk = 1:Nkpt
        for ik = 1:NkptTotal
            if umap[ikk] == mapping[ik]
                append!( list_ir_k, [kgrid[:,ik]] )
                break
            end
        end
    end

    kred = zeros(Float64,3,Nkpt)
    for ik = 1:Nkpt
        kred[1,ik] = list_ir_k[ik][1] / mesh[1]
        kred[2,ik] = list_ir_k[ik][2] / mesh[2]
        kred[3,ik] = list_ir_k[ik][3] / mesh[3]
    end
    
    # count for occurence of each unique mapping
    kcount = zeros(Int64,Nkpt)
    for ik = 1:Nkpt
        kcount[ik] = count( i -> ( i == umap[ik] ), mapping )
    end

    # calculate the weights
    wk = kcount[:]/sum(kcount)

    # need to calculate this here because PWGrid instance is not passed
    RecVecs = 2*pi*inv(atoms.LatVecs')

    return KPoints( Nkpt, kred, wk, RecVecs )

end


function gen_MonkhorstPack( mesh::Array{Int64,1} )
    ik = 0
    kpts = zeros(Float64,3,prod(mesh))
    for k = 1:mesh[3]
    for j = 1:mesh[2]
    for i = 1:mesh[1]
        ik = ik + 1
        kpts[1,ik] = (2*i - mesh[1] - 1)/(2*mesh[1])
        kpts[2,ik] = (2*j - mesh[2] - 1)/(2*mesh[2])
        kpts[3,ik] = (2*k - mesh[3] - 1)/(2*mesh[3])
    end
    end
    end
    return kpts
end

