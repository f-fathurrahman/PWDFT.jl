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

    G_to_R!( pw, ctmp )
    
    ortho_check( ctmp )
    
    ortho_gram_schmidt!( ctmp )
    
    ortho_check( ctmp )

    dVol = pw.CellVolume/Npoints
    ctmp[:] = ctmp[:]/sqrt(dVol)

    Rhoe = CuArrays.zeros( Float64, Npoints )

    for ist in 1:Nstates
        @views psi = ctmp[:,ist]
        Rhoe[:] = Rhoe[:] + real( conj(psi) .* psi )
    end
    println("integ Rhoe = ", sum(Rhoe)*dVol)


    println(1.0/Npoints)

    println("Pass here")
end

main()