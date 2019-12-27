using CuArrays
using FFTW

#CuArrays.allowscalar(false)

function simple_test()

    A = rand(ComplexF64,64,64,64)
    d_A = cu(A)

    for i in 1:2
        println("Using CPU")
        @time fft(A)
        println("Using GPU")
        @time fft(d_A)
    end
end
#simple_test()

function test_plan()

    Ns = (60,60,60)
    Npoints = prod(Ns)

    planfw = plan_fft( zeros(ComplexF64,Ns) )
    d_planfw = plan_fft( cu(zeros(ComplexF64,Ns)) )

    dat1 = rand(ComplexF64, Npoints)
    d_dat1 = cu(dat1)

    ft_dat1 = similar(dat1)
    d_ft_dat1 = similar(d_dat1)

    for i in 1:2
        @time ft_dat1 = planfw*reshape(dat1, Ns)
        @time d_ft_dat1 = d_planfw*reshape(d_dat1, Ns)
    end

    ft_dat1_gpu = collect(d_ft_dat1)

    println( "avg sum diff      = ", sum(ft_dat1_gpu - ft_dat1)/Npoints )
    println( "avg sum real diff = ", sum(abs.(real(ft_dat1_gpu) - real(ft_dat1)))/Npoints )
    println( "avg sum imag diff = ", sum(abs.(imag(ft_dat1_gpu) - imag(ft_dat1)))/Npoints )


    println("Pass here")
end
test_plan()
