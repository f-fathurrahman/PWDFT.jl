using FFTW
using PWDFT

function write_GVectors( f::IOStream, gvec::GVectors )
    write(f, gvec.Ng)
    write(f, gvec.G)
    write(f, gvec.G2)
    write(f, gvec.idx_g2r)
end

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

function write_KPoints( f::IOStream, kpoints::KPoints )
    write(f, kpoints.Nkpt)
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

    read!(f, k)
    read!(f, wk)
    read!(f, RecVecs)

    return KPoints(Nkpt, k, wk, RecVecs)
end

function test_kpoints()
    atoms = Atoms(xyz_string_frac="""
    2

    H  0.0   0.0   0.0
    H  0.25  0.25  0.25
    """, LatVecs=gen_lattice_fcc(5.0))
    kpoints = KPoints( atoms, [4,4,4], [0,0,0])

    println(kpoints)

    f = open("TEMP_kpoints.dat", "w")
    write_KPoints(f, kpoints)
    close(f)

    f = open("TEMP_kpoints.dat", "r")
    kpoints2 = read_KPoints(f)
    close(f)

    println(kpoints2)
end

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

function test01()
    r1 = rand(3)
    r2 = rand(2)

    println("r1 = ", r1)
    println("r2 = ", r2)
    
    f = open("TEMP_write.dat", "w")
    write(f, size(r1)[1])
    write(f, r1)
    write(f, size(r2)[1])
    write(f, r2)
    close(f)

    f = open("TEMP_write.dat", "r")
    #
    aN1 = Array{Int64}(undef,1)  # use array
    read!(f, aN1)
    N1 = aN1[1]
    println("N1 = ", N1)
    r1_read = Array{Float64}(undef,N1)
    read!(f, r1_read)
    println("r1_read = ", r1_read)
    #
    aN2 = Array{Int64}(undef,1)
    read!(f,aN2)
    N2 = aN2[1]
    println("N2 = ", N2)
    r2_read = Array{Float64}(undef,N2)
    read!(f, r2_read)
    println("r2_read = ", r2_read)

    close(f)
end

# Test write of Array{Array{Float64,2},1}
function test02()
    A = Array{Array{Float64,2},1}(undef,3)
    Ngw = [4, 5, 6]
    for i = 1:3
        A[i] = rand(Ngw[i],4)
    end
    println(A)
    
    f = open("TEMP_write.dat", "w")
    for i = 1:3
        write(f, A[i])
    end
    close(f)

    Aread = Array{Array{Float64,2},1}(undef,3)
    f = open("TEMP_write.dat", "r")
    # allocate memory
    for i = 1:3
        Aread[i] = Array{Float64}(undef,Ngw[i],4)
        read!(f, Aread[i])
    end
    close(f)

    for i = 1:3
        diff = sum(A[i] - Aread[i])
        println("i = ", i, " diff = ", diff)
    end

    println("test02 is finished")
end

function test_pw()

    atoms = Atoms(xyz_string_frac="""
    2

    H  0.0   0.0   0.0
    H  0.25  0.25  0.25
    """, LatVecs=gen_lattice_fcc(5.0))

    pw = PWGrid( 15.0, gen_lattice_fcc(5.0),
                 kpoints=KPoints( atoms, [4,4,4], [0,0,0]) )
    
    f = open("TEMP_pw.dat", "w")
    #write_GVectors(f, pw.gvec)
    #write_GVectorsW(f, pw.gvecw)
    write_PWGrid(f, pw)
    close(f)

    println(pw)

    f = open("TEMP_pw.dat", "r")
    #gvec = read_GVectors(f)
    #gvecw = read_GVectorsW(f)
    pw2 = read_PWGrid(f)
    close(f)

    println(pw2)
end


#test_kpoints()
test_pw()
#test01()
#test02()
