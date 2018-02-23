"""
Test orthonormalization of wavefunction
"""
function ortho_check( psi::Array{Complex128,2}; dVol=1.0 )
    Nstates = size(psi)[2]
    @printf("\nNorm check:\n")
    for ist = 1:Nstates
        c = dot( psi[:,ist], psi[:,ist] ) * dVol
        if c.im <= eps()
            @printf("State: #%5d: %18.10f\n", ist, c.re)
        else
            @printf("State: #%5d: (%18.10f,%18.10f)\n", ist, c.re, c.im)
        end
    end
    @printf("\nOrtho check w.r.t state #1:\n")
    for ist = 2:Nstates
        c = dot( psi[:,ist], psi[:,1] ) * dVol
        if c.im <= eps()
            @printf("State: #%5d: %18.10f\n", ist, c.re)
        else
            @printf("State: #%5d: (%18.10f,%18.10f)\n", ist, c.re, c.im)
        end
    end
    @printf("\n")
end
