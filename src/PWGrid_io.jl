function write_GVectors( f::IOStream, gvec::GVectors )
    write(f, gvec.Ng)
    write(f, gvec.G)
    write(f, gvec.G2)
    write(f, gvec.idx_g2r)
end

# GVectors --------------------------------------------------------------------

function read_GVectors( f::IOStream )
    tmpInt = Array{Int64}(undef,1)
    read!(f, tmpInt)
    Ng = tmpInt[1]

    G = Array{Float64}(undef,3,Ng)
    G2 = Array{Float64}(undef,Ng)
    idx_g2r = Array{Int64}(undef,Ng)

    read!(f, G)
    read!(f, G2)
    read!(f, idx_g2r)

    return GVectors(Ng, G, G2, idx_g2r)
end

# GVectorsW --------------------------------------------------------------------

function write_GVectorsW( f::IOStream, gvecw::GVectorsW )
    Nkpt = gvecw.kpoints.Nkpt
    write(f, Nkpt) # also write Nkpt, to preserve the order of writing
    write(f, gvecw.Ngwx)
    write(f, gvecw.Ngw)
    for ik = 1:Nkpt
        write(f, gvecw.idx_gw2g[ik])
        write(f, gvecw.idx_gw2r[ik])
    end
    write_KPoints( f, gvecw.kpoints )
end

function read_GVectorsW( f::IOStream )
    tmpInt = Array{Int64}(undef,1)
    read!(f, tmpInt)
    Nkpt = tmpInt[1]

    read!(f, tmpInt)
    Ngwx = tmpInt[1]

    Ngw = Array{Int64}(undef,Nkpt)
    read!(f, Ngw)

    idx_gw2g = Array{Array{Int64,1},1}(undef,Nkpt)
    idx_gw2r = Array{Array{Int64,1},1}(undef,Nkpt)

    for ik = 1:Nkpt
        idx_gw2g[ik] = Array{Int64}(undef,Ngw[ik])
        idx_gw2r[ik] = Array{Int64}(undef,Ngw[ik])
        read!(f, idx_gw2g[ik])
        read!(f, idx_gw2r[ik])
    end

    kpoints = read_KPoints(f)

    return GVectorsW(Ngwx, Ngw, idx_gw2g, idx_gw2r, kpoints)
end

# PWGrid --------------------------------------------------------------------

function write_PWGrid( f::IOStream, pw::PWGrid )
    write(f, pw.ecutwfc)
    write(f, pw.ecutrho)
    write(f, pw.Ns[1])
    write(f, pw.Ns[2])
    write(f, pw.Ns[3])
    write(f, pw.LatVecs)
    write(f, pw.RecVecs)
    write(f, pw.CellVolume)
    write(f, pw.r)
    write_GVectors(f, pw.gvec)
    write_GVectorsW(f, pw.gvecw)
    # planfw and planbw is not written, they will be recalculated from pw.Ns
end

function read_PWGrid(f)
    tmpInt = Array{Int64}(undef,1)
    tmpFlt = Array{Float64}(undef,1)
    
    read!(f, tmpFlt)
    ecutwfc = tmpFlt[1]

    read!(f, tmpFlt)
    ecutrho = tmpFlt[1]

    read!(f, tmpInt)
    Ns1 = tmpInt[1]

    read!(f, tmpInt)
    Ns2 = tmpInt[1]

    read!(f, tmpInt)
    Ns3 = tmpInt[1]

    Ns = (Ns1,Ns2,Ns3)

    LatVecs = Array{Float64}(undef,3,3)
    read!(f, LatVecs)

    RecVecs = Array{Float64}(undef,3,3)
    read!(f, RecVecs)

    read!(f, tmpFlt)
    CellVolume = tmpFlt[1]

    r = Array{Float64}(undef,3,prod(Ns))
    read!(f, r)

    gvec = read_GVectors(f)
    gvecw = read_GVectorsW(f)

    # Recalculate planfw and planbw
    planfw = plan_fft( zeros(Ns) )
    planbw = plan_ifft( zeros(Ns) )

    return PWGrid( ecutwfc, ecutrho, Ns, LatVecs, RecVecs, CellVolume, r, gvec, gvecw,
                   planfw, planbw )
    
end

# Overloaded println

import Base: println

"""
Display some information about `pw::PWGrid`. This function calls
`println(gvec::GVectors)` and `println(gvecw::GVectorsW)`.
"""
function println( pw::PWGrid; header=true )
    if header
        @printf("\n")
        @printf("                                     ------\n")
        @printf("                                     PWGrid\n")
        @printf("                                     ------\n")
        @printf("\n")
    end
    LatVecs = pw.LatVecs
    RecVecs = pw.RecVecs
    @printf("Direct lattice vectors:\n")
    @printf("\n")
    for i = 1:3
        @printf("%18.10f %18.10f %18.10f\n", LatVecs[i,1], LatVecs[i,2], LatVecs[i,3])
    end
    @printf("\n")
    @printf("Reciprocal lattice vectors:\n")
    @printf("\n")
    for i = 1:3
        @printf("%18.10f %18.10f %18.10f\n", RecVecs[i,1], RecVecs[i,2], RecVecs[i,3])
    end
    @printf("\n")
    @printf("Direct lattive volume = %18.10f bohr^3\n", pw.CellVolume )
    @printf("ecutwfc               = %18.10f Ha\n", pw.ecutwfc)
    @printf("ecutrho               = %18.10f Ha\n", pw.ecutrho)    
    @printf("Sampling points       = (%5d,%5d,%5d)\n", pw.Ns[1], pw.Ns[2], pw.Ns[3])
    #
    println( pw.gvec )
    println( pw.gvec, pw.gvecw )
end

"""
Display some information about `gvec::GVectors`.
"""
function println( gvec::GVectors )
    Ng = gvec.Ng
    G = gvec.G
    G2 = gvec.G2
    
    @printf("\n")
    @printf("                                    --------\n")
    @printf("                                    GVectors\n")
    @printf("                                    --------\n")
    @printf("\n")
    @printf("Ng = %12d\n", Ng)
    @printf("\n")
    for ig = 1:3
        @printf("%8d [%18.10f,%18.10f,%18.10f] : %18.10f\n", ig, G[1,ig], G[2,ig], G[3,ig], G2[ig])        
    end
    @printf(" ....... \n")
    for ig = Ng-3:Ng
        @printf("%8d [%18.10f.%18.10f,%18.10f] : %18.10f\n", ig, G[1,ig], G[2,ig], G[3,ig], G2[ig])
    end
    @printf("\n")
    @printf("Max G2 = %18.10f\n", maximum(G2))
end

"""
Display some information about `gvecw::GVectorsW`.
"""
function println( gvec::GVectors, gvecw::GVectorsW )
    G = gvec.G
    G2 = gvec.G2

    Ngwx = gvecw.Ngwx
    Ngw = gvecw.Ngw
    k = gvecw.kpoints.k
    Nkpt = gvecw.kpoints.Nkpt

    @printf("\n")
    @printf("                                    ---------\n")
    @printf("                                    GVectorsW\n")
    @printf("                                    ---------\n")
    @printf("\n")
    @printf("Ngwx = %12d\n", Ngwx)
        
    for ik = 1:Nkpt
        idx_gw2g = gvecw.idx_gw2g[ik]
        Gw = zeros(3,Ngw[ik])
        Gw2 = zeros(Ngw[ik])
        for igk = 1:Ngw[ik]
            ig = idx_gw2g[igk]
            Gw[:,igk] = G[:,ig] + k[:,ik]
            Gw2[igk] = Gw[1,igk]^2 + Gw[2,igk]^2 + Gw[3,igk]^2
        end
        @printf("Ngw = %8d, Max Gw2 = %18.10f\n", Ngw[ik], maximum(Gw2))
    end
end

