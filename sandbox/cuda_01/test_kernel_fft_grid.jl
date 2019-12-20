include("PWDFT_cuda.jl")

function main()

    pw = CuPWGrid( 15.0, gen_lattice_fcc(15.0) )

    Npoints = prod(pw.Ns)
    Nstates = 10
    Nspin = 1

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

    idx_cpu = collect(idx)
    ctmp_cpu[idx_cpu,:] = psi_cpu

    ctmp_gpu = collect(ctmp)

    ISTATES = 2
    for ip in 1:10
        @printf("%18.10f %18.10f\n", real(ctmp_cpu[ip,ISTATES]), real(ctmp_gpu[ip,ISTATES]))
    end

    println("Pass here")
end

main()