
# Overloaded show's

import Base: print

"""
Display some information about `pw::PWGrid`. This function calls
`println(gvec::GVectors)` and `println(gvecw::GVectorsW)`.
"""
function print( io::IO, pw::PWGrid; header=true )
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
    for i in 1:3
        @printf(io, "%18.10f %18.10f %18.10f\n", LatVecs[i,1], LatVecs[i,2], LatVecs[i,3])
    end
    @printf(io, "\n")
    @printf(io, "Reciprocal lattice vectors:\n")
    @printf(io, "\n")
    for i in 1:3
        @printf(io, "%18.10f %18.10f %18.10f\n", RecVecs[i,1], RecVecs[i,2], RecVecs[i,3])
    end
    @printf(io, "\n")
    @printf(io, "Direct lattive volume  = %18.10f bohr^3\n", pw.CellVolume )
    @printf(io, "ecutwfc                = %18.10f Ha\n", pw.ecutwfc)
    @printf(io, "ecutrho                = %18.10f Ha\n", pw.ecutrho)    
    @printf(io, "Sampling points        = (%5d,%5d,%5d)\n", pw.Ns[1], pw.Ns[2], pw.Ns[3])
    if pw.Nss != nothing
        @printf(io, "Sampling points smooth = (%5d,%5d,%5d)\n", pw.Nss[1], pw.Nss[2], pw.Nss[3])
    end
    #
    print( io, pw.gvec )
    if pw.using_dual_grid
        println("\nUsing dual grid")
        print( io, pw.gvecs )
        print( io, pw.gvecs, pw.gvecw )
    else
        print( io, pw.gvec, pw.gvecw )
    end
end
print( pw::PWGrid; header=true ) = print( stdout, pw, header=header )

"""
Display some information about `gvec::GVectors`.
"""
function print( io::IO, gvec::GVectors )
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
    for ig in 1:3
        @printf(io, "%8d [%18.10f,%18.10f,%18.10f] : %18.10f\n", ig, G[1,ig], G[2,ig], G[3,ig], G2[ig])        
    end
    @printf(io, " ....... \n")
    for ig in Ng-3:Ng
        @printf(io, "%8d [%18.10f.%18.10f,%18.10f] : %18.10f\n", ig, G[1,ig], G[2,ig], G[3,ig], G2[ig])
    end
    @printf(io, "\n")
    @printf(io, "Max G2 = %18.10f\n", maximum(G2))
end
print( gvec::GVectors ) = print( stdout, gvec )


"""
Display some information about `gvecw::GVectorsW`.
"""
function print( io::IO, gvec::GVectors, gvecw::GVectorsW )
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

    # FIXME: this loop does some calculations and allocates memory, remove it
    for ik in 1:Nkpt
        idx_gw2g = gvecw.idx_gw2g[ik]
        Gw = zeros(3,Ngw[ik])
        Gw2 = zeros(Ngw[ik])
        for igk in 1:Ngw[ik]
            ig = idx_gw2g[igk]
            @views Gw[:,igk] = G[:,ig] + k[:,ik]
            Gw2[igk] = Gw[1,igk]^2 + Gw[2,igk]^2 + Gw[3,igk]^2
        end
        @printf(io, "Ngw = %8d, Max Gw2 = %18.10f\n", Ngw[ik], maximum(Gw2))
        #@printf(io, "gvecw.Ngw[ik] = %8d (should be the same as Ngw above)\n", gvecw.Ngw[ik])
    end
end
print( gvec::GVectors, gvecw::GVectorsW ) = print( stdout, gvec, gvecw )

