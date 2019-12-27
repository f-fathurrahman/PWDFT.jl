using CUDAnative
using CuArrays
using Printf
using LinearAlgebra

using PWDFT
using PWDFT_cuda

function main()

    Nkpt = 4
    Ngw = [1000, 1100, 1200, 1200]
    Nstates = 5
    
    TYP = typeof( CuArray(rand(ComplexF64,1,1)) )

    psiks = Array{CuArray{ComplexF64,2},1}(undef,Nkpt)
    for ik in 1:Nkpt
        psiks[ik] = CuArray( rand(ComplexF64,Ngw[ik],Nstates) )
        ortho_gram_schmidt!( psiks[ik] )
    end

    ortho_check( psiks[1] )

    println(TYP)
    println(typeof(psiks[1]))

    println("Pass here")
end

main()