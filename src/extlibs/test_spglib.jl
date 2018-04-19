const SPGLIB_SO_PATH = "/home/efefer/WORKS/my_github_repos/PWDFT.jl/src/extlibs/spglib/libspglib.so"

function spg_find_primitive(lattice, position, types, num_atom, symprec)
    num_primitive_atom =
    ccall( (:my_spg_find_primitive,SPGLIB_SO_PATH), Cint,
           ( Ptr{Float64}, Ptr{Float64}, Ptr{Cint}, Cint, Float64 ),
           lattice, position, types, num_atom, symprec )
    return num_primitive_atom
end

function spg_find_primitive_v2(lattice, position, types, num_atom, symprec)
    println("Using spg_find_primitive_v2:")
    
    ctypes_cnv = Base.cconvert(Array{Int32,1},types)

    num_primitive_atom =
    ccall( (:spg_find_primitive,SPGLIB_SO_PATH), Int32,
           ( Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Int32, Float64 ),
           lattice, position, ctypes_cnv, num_atom, symprec )
    return num_primitive_atom
end

function test_call_2d_array()
    a = ones(Int32,3,3)
    println("a = ", a)
    ccall( (:test_2d_array, SPGLIB_SO_PATH), Void,
           (Ptr{Int32}, Int32), a, 3 )
end

#test_call_2d_array()

function test_BCC()
    lattice = 4.0*diagm(ones(3))

    Natoms = 2
    atpos = zeros(3,Natoms)
    atpos[:,1] = [0.0, 0.0, 0.0]
    atpos[:,2] = [0.5, 0.5, 0.5]

    types = zeros(Int64,2)
    types[1] = 1
    types[2] = 1

    println(typeof(types))
    symprec = 1e-5

    println(lattice)
    println(atpos)    

    println("Before num_primitive_atom")
    num_primitive_atom = spg_find_primitive_v2(lattice, atpos, types, Natoms, symprec)

    println("num_primitive_atom = ", num_primitive_atom)
end


function test_corrundum()
    lattice = zeros(3,3)

    lattice[1,:] = [4.8076344022756095, -2.4038172011378047, 0]
    lattice[2,:] = [0.0, 4.1635335244786962, 0.0]
    lattice[3,:] = [0.0, 0.0, 13.1172699198127543]

    Natoms = 30
    atpos = zeros(3,Natoms)
    atpos[:,1]  = [0.0000000000000000, 0.0000000000000000, 0.3521850942289043]
    atpos[:,2]  = [0.6666666666666643, 0.3333333333333357, 0.6855184275622400]
    atpos[:,3]  = [0.3333333333333357, 0.6666666666666643, 0.0188517608955686]
    atpos[:,4]  = [0.0000000000000000, 0.0000000000000000, 0.6478149057711028]
    atpos[:,5]  = [0.6666666666666643, 0.3333333333333357, 0.9811482391044314]
    atpos[:,6]  = [0.3333333333333357, 0.6666666666666643, 0.3144815724377600]
    atpos[:,7]  = [0.0000000000000000, 0.0000000000000000, 0.1478149057710957]
    atpos[:,8]  = [0.6666666666666643, 0.3333333333333357, 0.4811482391044314]
    atpos[:,9]  = [0.3333333333333357, 0.6666666666666643, 0.8144815724377600]
    atpos[:,10] = [0.0000000000000000, 0.0000000000000000, 0.8521850942288972]
    atpos[:,11] = [0.6666666666666643, 0.3333333333333357, 0.1855184275622400]
    atpos[:,12] = [0.3333333333333357, 0.6666666666666643, 0.5188517608955686]
    atpos[:,13] = [0.3061673906454899, 0.0000000000000000, 0.2500000000000000]
    atpos[:,14] = [0.9728340573121541, 0.3333333333333357, 0.5833333333333357]
    atpos[:,15] = [0.6395007239788255, 0.6666666666666643, 0.9166666666666643]
    atpos[:,16] = [0.6938326093545102, 0.0000000000000000, 0.7500000000000000]
    atpos[:,17] = [0.3604992760211744, 0.3333333333333357, 0.0833333333333357]
    atpos[:,18] = [0.0271659426878458, 0.6666666666666643, 0.4166666666666643]
    atpos[:,19] = [0.0000000000000000, 0.3061673906454899, 0.2500000000000000]
    atpos[:,20] = [0.6666666666666643, 0.6395007239788255, 0.5833333333333357]
    atpos[:,21] = [0.3333333333333357, 0.9728340573121541, 0.9166666666666643]
    atpos[:,22] = [0.0000000000000000, 0.6938326093545102, 0.7500000000000000]
    atpos[:,23] = [0.6666666666666643, 0.0271659426878458, 0.0833333333333357]
    atpos[:,24] = [0.3333333333333357, 0.3604992760211744, 0.4166666666666643]
    atpos[:,25] = [0.6938326093545102, 0.6938326093545102, 0.2500000000000000]
    atpos[:,26] = [0.3604992760211744, 0.0271659426878458, 0.5833333333333357]
    atpos[:,27] = [0.0271659426878458, 0.3604992760211744, 0.9166666666666643]
    atpos[:,28] = [0.3061673906454899, 0.3061673906454899, 0.7500000000000000]
    atpos[:,29] = [0.9728340573121541, 0.6395007239788255, 0.0833333333333357]
    atpos[:,30] = [0.6395007239788255, 0.9728340573121541, 0.4166666666666643]

    types = zeros(Cint,Natoms)
    symprec = 1e-5
    
    types[1:12]  = 1
    types[13:30] = 2

    println("Before find_primitive = ")
    for i = 1:3
        @printf("%f %f %f\n", lattice[i,1], lattice[i,2], lattice[i,3])
    end
    println("")    
    for ia = 1:Natoms
        @printf("%5d %18.10f %18.10f %18.10f\n", ia, atpos[1,ia], atpos[2,ia], atpos[3,ia])
    end    

    num_primitive_atom = spg_find_primitive_v2(lattice, atpos, types, Natoms, symprec)
    
    println("\nAfter find_primitive = ")
    println("")
    for i = 1:3
        @printf("%f %f %f\n", lattice[i,1], lattice[i,2], lattice[i,3])
    end
    println("")    
    for ia = 1:num_primitive_atom
        @printf("%5d %18.10f %18.10f %18.10f\n", ia, atpos[1,ia], atpos[2,ia], atpos[3,ia])
    end
    println("num_primitive_atom = ", num_primitive_atom)

end

#test_BCC()
test_corrundum()