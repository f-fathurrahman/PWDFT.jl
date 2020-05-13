function setup_guess_wavefunc!(
    Ham::HamiltonianGamma,
    psis::BlochWavefuncGamma,
    startingrhoe,
    skip_initial_diag
)
    
    Npoints = prod(Ham.pw.Ns)
    Nspin = Ham.electrons.Nspin
    Rhoe = zeros(Float64,Npoints,Nspin)
    
    if startingrhoe == :gaussian
        @assert Nspin == 1
        Rhoe[:,1] = guess_rhoe( Ham )
    else
        calc_rhoe!( Ham, psis, Rhoe )
    end
    #
    update!(Ham, Rhoe)
    #
    if !skip_initial_diag
        _ =
        diag_LOBPCG!( Ham, psis, verbose=false, verbose_last=false, NiterMax=10 )
    end
    return
end