import PWDFT: ortho_check

function ortho_check( psi::CuArray{ComplexF64,2}; dVol=1.0 )
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
