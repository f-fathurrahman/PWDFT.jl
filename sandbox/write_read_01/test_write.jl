using FFTW
using PWDFT

include("PWDFT_io.jl")

function test_kpoints()
    atoms = Atoms(xyz_string_frac="""
    2

    H  0.0   0.0   0.0
    H  0.25  0.25  0.25
    """, LatVecs=gen_lattice_fcc(5.0))
    kpoints = KPoints(atoms, [4,4,4], [0,0,0])

    println(kpoints)

    f = open("TEMP_kpoints.dat", "w")
    write_KPoints(f, kpoints)
    close(f)

    f = open("TEMP_kpoints.dat", "r")
    kpoints2 = read_KPoints(f)
    close(f)

    println(kpoints2)
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


# writing Array{String}
function test03()

    atoms = Atoms(xyz_string_frac="""
    2

    H  0.0   0.0   0.0
    H  0.25  0.25  0.25
    """, LatVecs=gen_lattice_fcc(5.0))

    f = open("TEMP_string.dat", "w")
    
    write(f, atoms.positions)
    for ia = 1:atoms.Natoms
        write(f, codeunits(atoms.atsymbs[ia]))
    end
    write(f, atoms.Zvals)

    close(f)

    f = open("TEMP_string.dat", "r")
    positions = Array{Float64,2}(undef,3,atoms.Natoms)
    read!(f, positions)

    tmpStr = Array{UInt8}(undef,1)
    read!(f, tmpStr)

    Zvals = Array{Float64}(undef,atoms.Nspecies)
    read!(f, Zvals)

    close(f)

    println("positions = ", positions)
    println("tmpStr = ", tmpStr)
    println("Zvals = ", Zvals)

end
test03()

function test_pw()

    atoms = Atoms(xyz_string_frac="""
    2

    H  0.0   0.0   0.0
    H  0.25  0.25  0.25
    """, LatVecs=gen_lattice_fcc(5.0))

    pw = PWGrid( 15.0, gen_lattice_fcc(5.0),
                 kpoints=KPoints( atoms, [4,4,4], [0,0,0]) )
    
    f = open("TEMP_pw.dat", "w")
    write_PWGrid(f, pw)
    close(f)

    println(pw)

    f = open("TEMP_pw.dat", "r")
    pw2 = read_PWGrid(f)
    close(f)

    println(pw2)
end


#test_kpoints()
#test_pw()
#test01()
#test02()


function test_Atoms()
    atoms = Atoms(xyz_string_frac="""
    2

    H  0.0   0.0   0.0
    H  0.25  0.25  0.25
    """, LatVecs=gen_lattice_fcc(5.0))

    f = open("TEMP_atoms.dat", "w")
    write_Atoms(f, atoms)
    close(f)

    println(atoms)

    f = open("TEMP_pw.dat", "r")
    atoms2 = read_Atoms(f)
    close(f)

    println(atoms2)
end
#test_Atoms()
