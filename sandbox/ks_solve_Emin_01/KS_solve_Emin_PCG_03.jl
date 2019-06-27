function KS_solve_Emin_PCG_03!(
    Ham::Hamiltonian;
    startingwfc=:random,
    savewfc=false,
    startingrhoe=:gaussian,
    skip_initial_diag=false,
    α_t=3e-5,
    NiterMax=200,
    verbose=true,
    print_final_ebands=false,
    print_final_energies=true,
    i_cg_beta=2,
    etot_conv_thr=1e-6,
    kT=0.001
)


    Haux = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    eta = zeros(Float64,Nstates,Nkspin)
    for ikspin = 1:Nkspin
        Haux[ikspin] = rand(ComplexF64, Nstates, Nstates)
        Haux[ikspin] = 0.5*(Haux[ikspin] + Haux[ikspin]')
        eta[:,ikspin] = eigvals(Haux[ikspin])
    end

    Ham.electrons.Focc[:,:], E_fermi = calc_Focc( Nelectrons, wk, kT, eta, Nspin )
    println("E_fermi = ", E_fermi)
    for ist = 1:Nstates
        println(eta[ist,1], " ", Ham.electrons.Focc[ist,1])
    end

    Hsub = Array{Array{ComplexF64,2},1}(undef,Nkspin)

    g_Haux = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    gt_Haux = Array{Array{ComplexF64,2},1}(undef,Nkspin)    
    Kg_Haux = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    g_Haux_old = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    Kg_Haux_old = Array{Array{ComplexF64,2},1}(undef,Nkspin)    
    Haux_c = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    d_Haux = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    d_Haux_old = Array{Array{ComplexF64,2},1}(undef,Nkspin)

    for ikspin = 1:Nkspin

        Hsub[ikspin] = zeros(ComplexF64,Nstates,Nstates)

        g_Haux[ikspin] = zeros(ComplexF64,Nstates,Nstates)
        gt_Haux[ikspin] = zeros(ComplexF64,Nstates,Nstates)        
        Kg_Haux[ikspin] = zeros(ComplexF64,Nstates,Nstates)
        g_Haux_old[ikspin] = zeros(ComplexF64,Nstates,Nstates)
        Kg_Haux_old[ikspin] = zeros(ComplexF64,Nstates,Nstates)        
        Haux_c[ikspin] = zeros(ComplexF64,Nstates,Nstates)
        d_Haux[ikspin] = zeros(ComplexF64,Nstates,Nstates)
        d_Haux_old[ikspin] = zeros(ComplexF64,Nstates,Nstates)        
    end


    


    return

end