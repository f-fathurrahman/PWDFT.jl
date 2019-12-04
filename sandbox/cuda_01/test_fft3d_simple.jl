using CuArrays
using FFTW
using Printf

function main( ; NN=(2,2,2) )

    A = rand(ComplexF64, NN[1], NN[2], NN[3])

    d_A = cu(A)
    d_A_G = fft(d_A)
    A_G = collect(d_A_G)

    A_G_ref = fft(A)

    N = length(A)
    if N < 10
        for i in 1:N
            @printf("%18.10f %18.10f %18.10f\n", abs(A_G[i]), abs(A_G_ref[i]), abs(A_G[i]-A_G_ref[i]))
        end
    end
    println("diff = ", sum(A_G - A_G_ref)/N)
end

main(NN=(10,10,10))
