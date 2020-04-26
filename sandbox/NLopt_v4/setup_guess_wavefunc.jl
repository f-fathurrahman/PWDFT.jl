function setup_guess_wavefunc!( Ham, psiks; startingrhoe=:gaussian, skip_initial_diag=false )
    
    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Rhoe = zeros(Float64,Npoints,Nspin)
    
    if startingrhoe == :gaussian
        @assert Nspin == 1
        Rhoe[:,1] = guess_rhoe( Ham )
    else
        calc_rhoe!( Ham, psiks, Rhoe )
    end
    #
    update!(Ham, Rhoe)
    #
    if !skip_initial_diag
        _ =
        diag_LOBPCG!( Ham, psiks, verbose=false, verbose_last=false, NiterMax=10 )
    end
    return
end