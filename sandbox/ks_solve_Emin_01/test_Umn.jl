using PWDFT
using LinearAlgebra

function main()
    Nstates = 4
    Nbasis = 100
    
    psi = rand(ComplexF64,Nbasis,Nstates)
    ortho_check(psi)

    U = inv(sqrt(psi'*psi))

    psi_ortho = psi*U
    ortho_check(psi_ortho)

    eta = rand(ComplexF64,Nstates,Nstates)
    print_matrix(U'*eta*U)
end

main()
