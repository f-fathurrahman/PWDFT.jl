"""
Test orthonormalization of wavefunction
"""
function ortho_check( psi::Array{ComplexF64,2}; dVol=1.0 )
    Nstates = size(psi)[2]
    @printf("\nNorm check:\n")
    for ist = 1:Nstates
        c = dot( psi[:,ist], psi[:,ist] ) * dVol
        @printf("State: #%5d: (%18.10f,%18.10f)\n", ist, c.re, c.im)
    end
    @printf("\nOrtho check w.r.t state #1:\n")
    for ist = 2:Nstates
        c = dot( psi[:,ist], psi[:,1] ) * dVol
        @printf("State: #%5d: (%18.10f,%18.10f)\n", ist, c.re, c.im)
    end
    @printf("\n")
end


function ortho_check( Ham::Hamiltonian, psi::Array{ComplexF64,2}; dVol=1.0 )
    Nstates = size(psi,2)
    if Ham.need_overlap
        Spsi = op_S(Ham, psi)
    else
        Spsi = psi
    end
    @printf("\nNorm check:\n")
    for ist in 1:Nstates
        c = dot( psi[:,ist], Spsi[:,ist] ) * dVol
        @printf("State: #%5d: (%18.10f,%18.10f)\n", ist, c.re, c.im)
    end
    @printf("\nOrtho check w.r.t state #1:\n")
    for ist in 1:Nstates
        c = dot( psi[:,ist], Spsi[:,1] ) * dVol
        @printf("State: #%5d: (%18.10f,%18.10f)\n", ist, c.re, c.im)
    end
    @printf("\n")
    return
end