using Printf
using PWDFT

function test_no_spin()
    
    Nstates = 5
    Nkpt = 1
    Nspin = 1
    Nkspin = Nkpt*Nspin
    kT = 0.01

    Focc = zeros(Nstates,Nkspin)
    Focc[:,1] = [2.0, 2.0, 1.8, 0.2, 0.0]
    
    evals = zeros(Nstates,Nkspin)
    evals[:,1] = [-2.0, -1.0, -0.5, 0.5, 1.0]
    
    E_fermi = 0.1
    wk = ones(Nkpt)/Nkpt

    mTS = calc_entropy( Focc, wk, kT, evals, E_fermi )
    @printf("mTS = %18.10f\n", mTS)

    mTS = calc_entropy_v2( wk, kT, evals, E_fermi, Nspin )
    @printf("mTS v2 = %18.10f\n", mTS)

end

function test_with_spin()
    
    Nstates = 5
    Nkpt = 1
    Nspin = 2
    Nkspin = Nkpt*Nspin
    kT = 0.01

    Focc = zeros(Nstates,Nkspin)
    Focc[:,1] = [1.0, 1.0, 0.9, 0.1, 0.0]
    Focc[:,2] = [1.0, 1.0, 0.9, 0.1, 0.0]
    
    evals = zeros(Nstates,Nkspin)
    evals[:,1] = [-2.0, -1.0, -0.5, 0.5, 1.0]
    evals[:,2] = [-2.0, -1.0, -0.5, 0.5, 1.0]
    
    E_fermi = 0.1
    wk = ones(Nkpt)/Nkpt

    mTS = calc_entropy( Focc, wk, kT, evals, E_fermi, Nspin=Nspin )
    @printf("mTS = %18.10f\n", mTS)
end


test_no_spin()
test_with_spin()

