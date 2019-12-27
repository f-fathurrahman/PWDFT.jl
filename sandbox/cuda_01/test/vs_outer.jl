using CuArrays
using Printf

function main(; Nbasis=10, Nstates=5)
    psi = rand( ComplexF64, Nbasis, Nstates )
    d_psi = CuArray(psi)

    pp = psi' * psi

    @printf("Using CPU:                      ")
    @time pp = psi' * psi

    d_pp = d_psi' * d_psi

    @printf("Using GPU:                      ")
    @time d_pp = d_psi' * d_psi

    @printf("Using GPU + copy back the result")
    @time begin
        d_pp = d_psi' * d_psi
        pp_v2 = collect(d_pp)
    end

    if Nstates <= 5
        @printf("diff = %18.10f\n", sum(abs.(pp - pp_v2))/length(pp))
    end
end

main(Nbasis=10000, Nstates=50)
