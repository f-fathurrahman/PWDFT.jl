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

function write_KPoints( f::IOStream, kpoints::KPoints )
    write(f, kpoints.Nkpt)
    write(f, kpoints.mesh[1])
    write(f, kpoints.mesh[2])
    write(f, kpoints.mesh[3])
    write(f, kpoints.k)
    write(f, kpoints.wk)
    write(f, kpoints.RecVecs)
end

function read_KPoints( f::IOStream )
    
    tmpInt = Array{Int64}(undef,1)
    
    read!(f, tmpInt)
    Nkpt = tmpInt[1]
    
    k = Array{Float64}(undef,3,Nkpt)
    wk = Array{Float64}(undef,Nkpt)
    RecVecs = Array{Float64}(undef,3,3)

    read!(f, tmpInt)
    mesh1 = tmpInt[1]
    
    read!(f, tmpInt)
    mesh2 = tmpInt[1]
    
    read!(f, tmpInt)
    mesh3 = tmpInt[1]
    
    mesh = (mesh1, mesh2, mesh3)

    read!(f, k)
    read!(f, wk)
    read!(f, RecVecs)

    return KPoints(Nkpt, mesh, k, wk, RecVecs)
end