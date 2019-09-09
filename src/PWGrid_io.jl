
# Overloaded println

import Base: show

"""
Display some information about `pw::PWGrid`. This function calls
`println(gvec::GVectors)` and `println(gvecw::GVectorsW)`.
"""
function show( io::IO, pw::PWGrid; header=true )
    if header
        @printf("\n")
        @printf("                                     ------\n")
        @printf("                                     PWGrid\n")
        @printf("                                     ------\n")
        @printf("\n")
    end
    LatVecs = pw.LatVecs
    RecVecs = pw.RecVecs
    @printf(io, "Direct lattice vectors:\n")
    @printf(io, "\n")
    for i = 1:3
        @printf(io, "%18.10f %18.10f %18.10f\n", LatVecs[i,1], LatVecs[i,2], LatVecs[i,3])
    end
    @printf(io, "\n")
    @printf(io, "Reciprocal lattice vectors:\n")
    @printf(io, "\n")
    for i = 1:3
        @printf(io, "%18.10f %18.10f %18.10f\n", RecVecs[i,1], RecVecs[i,2], RecVecs[i,3])
    end
    @printf(io, "\n")
    @printf(io, "Direct lattive volume = %18.10f bohr^3\n", pw.CellVolume )
    @printf(io, "ecutwfc               = %18.10f Ha\n", pw.ecutwfc)
    @printf(io, "ecutrho               = %18.10f Ha\n", pw.ecutrho)    
    @printf(io, "Sampling points       = (%5d,%5d,%5d)\n", pw.Ns[1], pw.Ns[2], pw.Ns[3])
    #
    show( io, pw.gvec )
    show( io, pw.gvec, pw.gvecw )
end
show( pw::PWGrid ) = show( stdout, pw )

"""
Display some information about `gvec::GVectors`.
"""
function show( io::IO, gvec::GVectors )
    Ng = gvec.Ng
    G = gvec.G
    G2 = gvec.G2
    
    @printf(io, "\n")
    @printf(io, "                                    --------\n")
    @printf(io, "                                    GVectors\n")
    @printf(io, "                                    --------\n")
    @printf(io, "\n")
    @printf(io, "Ng = %12d\n", Ng)
    @printf(io, "\n")
    for ig = 1:3
        @printf(io, "%8d [%18.10f,%18.10f,%18.10f] : %18.10f\n", ig, G[1,ig], G[2,ig], G[3,ig], G2[ig])        
    end
    @printf(io, " ....... \n")
    for ig = Ng-3:Ng
        @printf(io, "%8d [%18.10f.%18.10f,%18.10f] : %18.10f\n", ig, G[1,ig], G[2,ig], G[3,ig], G2[ig])
    end
    @printf(io, "\n")
    @printf(io, "Max G2 = %18.10f\n", maximum(G2))
end
show( gvec::GVectors ) = show( stdout, gvec )


"""
Display some information about `gvecw::GVectorsW`.
"""
function show( io::IO, gvec::GVectors, gvecw::GVectorsW )
    G = gvec.G
    G2 = gvec.G2

    Ngwx = gvecw.Ngwx
    Ngw = gvecw.Ngw
    k = gvecw.kpoints.k
    Nkpt = gvecw.kpoints.Nkpt

    @printf(io, "\n")
    @printf(io, "                                    ---------\n")
    @printf(io, "                                    GVectorsW\n")
    @printf(io, "                                    ---------\n")
    @printf(io, "\n")
    @printf(io, "Ngwx = %12d\n", Ngwx)
        
    for ik = 1:Nkpt
        idx_gw2g = gvecw.idx_gw2g[ik]
        Gw = zeros(3,Ngw[ik])
        Gw2 = zeros(Ngw[ik])
        for igk = 1:Ngw[ik]
            ig = idx_gw2g[igk]
            Gw[:,igk] = G[:,ig] + k[:,ik]
            Gw2[igk] = Gw[1,igk]^2 + Gw[2,igk]^2 + Gw[3,igk]^2
        end
        @printf(io, "Ngw = %8d, Max Gw2 = %18.10f\n", Ngw[ik], maximum(Gw2))
    end
end
show( gvec::GVectors, gvecw::GVectorsW ) = show( stdout, gvec, gvecw )

