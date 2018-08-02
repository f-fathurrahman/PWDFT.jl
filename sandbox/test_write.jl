function test_main()
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

test_main()

