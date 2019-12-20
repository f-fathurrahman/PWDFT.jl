include("PWDFT_cuda.jl")

function main()

    Nbasis  = 10
    Nstates = 4

    psi = CuArrays.rand( ComplexF64, Nbasis, Nstates )
    psi_orig = psi[:,:]

    ortho_check( psi )
    ortho_gram_schmidt!( psi )
    ortho_check( psi )

    psi_cpu = collect(psi_orig)
    
    ortho_check( psi_cpu )
    #ortho_sqrt!( psi_cpu )
    ortho_gram_schmidt!( psi_cpu )
    ortho_check( psi_cpu )

    psi_gpu = collect( psi )
    ist = 1
    println("Some psi for ist = ", ist)
    for i = 1:5
        @printf("%18.10f %18.10f\n", real(psi_cpu[i,ist]), real(psi_gpu[i,ist]))
    end

end

main()
