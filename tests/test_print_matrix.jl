using PWDFT: print_matrix

function test_main()
    A = rand(Float64,3,3)
    println("\nMatrix A:")
    print_matrix(A)
    #
    B = rand(Float64,3,3) + im*rand(Float64,3,3)
    println("\nMatrix B:")
    print_matrix(B)
end

test_main()
