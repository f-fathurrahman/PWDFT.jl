using PWDFT

function test_main()
    srand(2345)
    psi = rand(1000,4) + im*rand(1000,4)
    #
    println("\nBefore ortho_gram_schmidt:")
    ortho_check(psi)
    #
    psi = ortho_gram_schmidt(psi)
    println("\nAfter ortho_gram_schmidt:")
    ortho_check(psi)
end


test_main()
