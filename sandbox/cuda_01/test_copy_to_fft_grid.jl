include("PWDFT_cuda.jl")
include("cu_op_V_loc.jl")

function main()

    pw = CuPWGrid( 15.0, gen_lattice_fcc(15.0) )

    Npoints = prod(pw.Ns)
    Nstates = 10
    Nspin = 1

    ctmp = CuArrays.fill( 0.0 + im*0.0, (Npoints,Nstates) )
    
    psiks = rand_CuBlochWavefunc( pw, Nstates, Nspin )

    ik = 1
    Nthreads = 256
    Ngw = pw.gvecw.Ngw
    Nblocks = ceil(Int64, Ngw[ik]/Nthreads)
    
    println("Nblocks = ", Nblocks)

    idx = pw.gvecw.idx_gw2r[ik]
    for ist in 1:Nstates
        @views psi = psiks[ik][:,ist]
        @views cc = ctmp[:,ist]
        println( dot(psi,psi) )
        @cuda threads=Nthreads blocks=Nblocks kernel_copy_to_fft_grid_gw2r(idx, psi, cc )
    end

    G_to_R!( pw, ctmp )
    
    ortho_check( ctmp )
    
    ortho_gram_schmidt!( ctmp )
    
    ortho_check( ctmp )

    dVol = pw.CellVolume/Npoints
    ctmp[:] = ctmp[:]/sqrt(dVol)

    Rhoe = CuArrays.fill( 0.0, Npoints )

    for ist in 1:Nstates
        @views psi = ctmp[:,ist]
        Rhoe[:] = Rhoe[:] + real( conj(psi) .* psi )
    end
    println("integ Rhoe = ", sum(Rhoe)*dVol)


    println(1.0/Npoints)

    println("Pass here")
end

main()