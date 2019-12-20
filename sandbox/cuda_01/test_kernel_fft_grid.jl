using Test

include("PWDFT_cuda.jl")

function main()

    ecutwfc = 15.0
    LatVecs = gen_lattice_fcc(10.0)

    pw = CuPWGrid( ecutwfc, LatVecs )

    pw_cpu = PWGrid( ecutwfc, LatVecs )

    Npoints = prod(pw.Ns)
    Nstates = 10
    Nspin = 1
    CellVolume = pw.CellVolume
    dVol = CellVolume/Npoints

    ctmp = CuArrays.zeros( ComplexF64, Npoints, Nstates )
    
    psiks = rand_CuBlochWavefunc( pw, Nstates, Nspin )

    ik = 1
    Nthreads = 256
    Ngw = pw.gvecw.Ngw
    Nblocks = ceil(Int64, Ngw[ik]/Nthreads)
    
    println("Nblocks = ", Nblocks)

    idx = pw.gvecw.idx_gw2r[ik]
    for ist in 1:Nstates
        @cuda threads=Nthreads blocks=Nblocks kernel_copy_to_fft_grid_gw2r!( ist, idx, psiks[ik], ctmp )
    end

    ctmp_cpu = zeros( ComplexF64, Npoints, Nstates )
    psi_cpu = collect(psiks[1])

    idx_cpu = pw_cpu.gvecw.idx_gw2r[ik]
    ctmp_cpu[idx_cpu,:] = psi_cpu[:,:]

    ctmp_gpu = collect(ctmp)

    ISTATES = 2
    println("Real part:")
    for ip in 1:10
        @printf("%18.10f %18.10f\n", real(ctmp_cpu[ip,ISTATES]), real(ctmp_gpu[ip,ISTATES]))
    end
    println("Imaginary:")
    for ip in 1:10
        @printf("%18.10f %18.10f\n", imag(ctmp_cpu[ip,ISTATES]), imag(ctmp_gpu[ip,ISTATES]))
    end

    #
    # GPU
    #
    # Bring to real space
    G_to_R!( pw, ctmp )
    ctmp[:] = sqrt(Npoints/CellVolume)*sqrt(Npoints)*ctmp[:]
    # Calculate Rhoe
    Rhoe = CuArrays.zeros(Float64,Npoints,Nspin)
    for ist = 1:Nstates
        Rhoe[:,1] = Rhoe[:,1] + real( conj(ctmp[:,ist]) .* ctmp[:,ist] )
    end
    println("integ Rhoe GPU (directly) = ", sum(Rhoe)*dVol)


    #
    # CPU
    #
    # Bring to real space
    G_to_R!( pw_cpu, ctmp_cpu )
    ctmp_cpu[:] = sqrt(Npoints/CellVolume)*sqrt(Npoints)*ctmp_cpu[:]
    # Calculate Rhoe
    Rhoe_cpu = zeros(Float64,Npoints,Nspin)
    for ist = 1:Nstates
        Rhoe_cpu[:,1] = Rhoe_cpu[:,1] + real( conj(ctmp_cpu[:,ist]) .* ctmp_cpu[:,ist] )
    end
    println("integ Rhoe CPU = ", sum(Rhoe_cpu)*dVol)


    # Compare
    Rhoe_gpu = collect(Rhoe)
    println("Some Rhoe")
    for ip in 1:10
        @printf("%18.10f %18.10f\n", Rhoe_cpu[ip,1], Rhoe_gpu[ip,1])
    end


    println("Pass here")
end

main()