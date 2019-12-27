using Printf

using CuArrays
using PWDFT
using PWDFT_cuda

function main()

    L = 20.0
    ecutwfc = 15.0

    pw = CuPWGrid( ecutwfc, gen_lattice_sc(L) )

    Npoints = prod(pw.Ns)
    Rhoe = CuArrays.rand( Float64, Npoints )

    phiG = Poisson_solve( pw, Rhoe )

    pw_ = PWGrid( ecutwfc, gen_lattice_sc(L) )
    Rhoe_ = collect( Rhoe )
    phiG_cpu = Poisson_solve( pw_, Rhoe_ )

    phiG_gpu = collect( phiG )

    println("diff = ", sum(pw_.gvec.G2 - collect(pw.gvec.G2)))
    println("diff = ", sum(Rhoe_ - collect(Rhoe)))
    
    println("diff = ", sum(phiG_cpu - phiG_gpu)/Npoints)

    for i in 1:5
        @printf("[%18.10f %18.10f], [%18.10f %18.10f]\n",
            real(phiG_cpu[i]), imag(phiG_cpu[i]),
            real(phiG_gpu[i]), imag(phiG_gpu[i]) )
    end

    for i in (Npoints-5):Npoints
        @printf("[%18.10f %18.10f], [%18.10f %18.10f]\n",
            real(phiG_cpu[i]), imag(phiG_cpu[i]),
            real(phiG_gpu[i]), imag(phiG_gpu[i]) )
    end

    println(typeof(phiG_gpu))
    println(typeof(phiG_cpu))
    println(typeof(phiG))

    println("Pass here")
end


main()
