using PWDFT
using PWDFT.PW

mutable struct PWHamiltonian
    pw::PWGrid
    V_ps_loc::Array{Float64,1}
    V_xc::Array{Float64,1}
    V_Hartree::Array{Float64,1}
    Rhoe::Array{Float64,1}
    Focc::Array{Float64,1}
end

function init_PWHamiltonian( pw, atoms::Atoms )
end

function op_H(H::PWHamiltonian, psi::Array{Complex128,2})
end

function test_main()

end