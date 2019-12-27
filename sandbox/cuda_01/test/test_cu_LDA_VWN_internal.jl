using Test
using Random
using Printf

using CuArrays
using PWDFT
using PWDFT_cuda

function main_nospin()
    Npoints = 1000
    Rhoe = abs.( CuArrays.rand(Float64, Npoints, 1) )

    epsxc = calc_epsxc_VWN( XCCalculator(), Rhoe )
    Vxc = calc_Vxc_VWN( XCCalculator(), Rhoe )

    epsxc_cpu = calc_epsxc_VWN( XCCalculator(), collect(Rhoe) )
    Vxc_cpu = calc_Vxc_VWN( XCCalculator(), collect(Rhoe) )

    epsxc_gpu = collect( epsxc )
    Vxc_gpu = collect( Vxc )

    println("Some epsxc")
    for i in 1:5
        @printf("%18.10f %18.10f\n", epsxc_cpu[i], epsxc_gpu[i])
    end

    println("Some Vxc")
    for i in 1:5
        @printf("%18.10f %18.10f\n", Vxc_cpu[i], Vxc_gpu[i])
    end

    @test epsxc_cpu ≈ epsxc_gpu
    @test Vxc_cpu ≈ Vxc_gpu

    println("Pass here")
end

function main_spin()
    Npoints = 1000
    Rhoe = abs.( CuArrays.rand(Float64, Npoints, 2) )

    epsxc = calc_epsxc_VWN( XCCalculator(), Rhoe )
    Vxc = calc_Vxc_VWN( XCCalculator(), Rhoe )

    epsxc_cpu = calc_epsxc_VWN( XCCalculator(), collect(Rhoe) )
    Vxc_cpu = calc_Vxc_VWN( XCCalculator(), collect(Rhoe) )

    epsxc_gpu = collect( epsxc )
    Vxc_gpu = collect( Vxc )

    println("Some epsxc")
    for i in 1:5
        @printf("%18.10f %18.10f\n", epsxc_cpu[i], epsxc_gpu[i])
    end

    println("Some Vxc spin up")
    for i in 1:5
        @printf("%18.10f %18.10f\n", Vxc_cpu[i,1], Vxc_gpu[i,1])
    end

    println("Some Vxc spin down")
    for i in 1:5
        @printf("%18.10f %18.10f\n", Vxc_cpu[i,2], Vxc_gpu[i,2])
    end

    @test epsxc_cpu ≈ epsxc_gpu
    @test Vxc_cpu ≈ Vxc_gpu

    println("Pass here")
end


#main_nospin()
main_spin()