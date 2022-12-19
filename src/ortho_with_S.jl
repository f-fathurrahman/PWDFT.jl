function ortho_sqrt_with_S!( Ham::Hamiltonian, psi::Array{ComplexF64,2} )
    O = psi' * op_S(Ham, psi)
    Udagger = inv(sqrt(O))
    psi[:,:] = psi*Udagger
    return
end


function ortho_check_with_S( Ham::Hamiltonian, psi::Array{ComplexF64,2}; dVol=1.0 )
    Nstates = size(psi,2)
    Spsi = op_S(Ham, psi)
    @printf("\nNorm check:\n")
    for ist = 1:Nstates
        c = dot( psi[:,ist], Spsi[:,ist] ) * dVol
        @printf("State: #%5d: (%18.10f,%18.10f)\n", ist, c.re, c.im)
    end
    @printf("\nOrtho check w.r.t state #1:\n")
    for ist = 1:Nstates
        c = dot( psi[:,ist], Spsi[:,2] ) * dVol
        @printf("State: #%5d: (%18.10f,%18.10f)\n", ist, c.re, c.im)
    end
    @printf("\n")
    return
end