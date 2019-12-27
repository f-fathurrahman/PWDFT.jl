using CuArrays
using CuArrays.CUSOLVER
using Printf
using LinearAlgebra

using PWDFT
using PWDFT_cuda

function main( ; Nbasis=1000, Nstates=10, Nkspin=5 )

    println()

    psiks = Array{CuArray{ComplexF64,2},1}(undef, Nkspin)

    @printf("Allocating array: ")
    @time begin
    for i in 1:Nkspin
        psiks[i] = CuArray( rand(ComplexF64,Nbasis,Nstates) )
    end
    end

    @printf("Operations with array: ")
    @time begin
    for i in 1:Nkspin
        psiks[i] = 2*psiks[i]
    end
    end

end

function main_noGPU( ; Nbasis=1000, Nstates=10, Nkspin=5 )

    println()

    psiks = Array{Array{ComplexF64,2},1}(undef, Nkspin)

    @printf("Allocating array: ")
    @time begin
    for i in 1:Nkspin
        psiks[i] = rand(ComplexF64, Nbasis,Nstates)
    end
    end

    @printf("Operations with array: ")
    @time begin
    for i in 1:Nkspin
        psiks[i] = 2*psiks[i]
    end
    end

end


function driver()

    Nbasis = 5000
    Nstates = 200
    Nkspin = 20

    println("No GPU")
    main_noGPU(Nbasis=Nbasis, Nstates=Nstates, Nkspin=Nkspin)
    main_noGPU(Nbasis=Nbasis, Nstates=Nstates, Nkspin=Nkspin)

    println()
    println("With GPU:")
    main(Nbasis=Nbasis, Nstates=Nstates, Nkspin=Nkspin)
    main(Nbasis=Nbasis, Nstates=Nstates, Nkspin=Nkspin)
end

driver()
