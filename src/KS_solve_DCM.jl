function KS_solve_DCM!( Ham::PWHamiltonian;
                        NiterMax = 100 )


	pw = Ham.pw
    Ngwx = pw.gvecw.Ngwx
    Ns = pw.Ns
    Npoints = prod(Ns)
    ΔV = pw.Ω/Npoints
    Focc = Ham.electrons.Focc
    Nstates = Ham.electrons.Nstates

    #
    # Random guess of wave function
    #
    srand(1234)
    psi = rand(Complex128,Ngwx,Nstates)
    psi = ortho_gram_schmidt(psi)
    #
    # Calculated electron density from this wave function and update Hamiltonian
    #
    rhoe = calc_rhoe( pw, Focc, psi )
    update!(Ham, rhoe)

    # Starting eigenvalues and psi
    λ, psi = diag_lobpcg( Ham, psi, verbose_last=false, maxit=10 )

    Energies = calc_energies( Ham, psi )
    Ham.energies = Energies
    Etot_old = Energies.Total

    # subspace
    Y = zeros( Complex128, Ngwx, 3*Nstates )
    R = zeros( Complex128, Ngwx, Nstates )
    P = zeros( Complex128, Ngwx, Nstates )

    set1 = 1:Nstates
    set2 = Nstates+1:2*Nstates
    set3 = 2*Nstates+1:3*Nstates
    set4 = Nstates+1:3*Nstates

    MaxInnerSCF = 3

    for iter = 1:NiterMax

        Hpsi = op_H( Ham, psi )
        #
        T = psi' * Hpsi
        T = 0.5*( T + T' )
        #
        R = Hpsi - psi*T
        R = Kprec( pw, R )
        #
        Y[:,set1] = psi
        Y[:,set2] = R
        if iter > 1
            Y[:,set3] = P
        end
        #
        # Project kinetic and ionic potential
        #
        KY = op_K( Ham, Y ) + op_V_Ps_loc( Ham, Y )
        T = real( Y'*KY )
        B = real( Y'*Y )
        B = 0.5*( B + B' )

        if iter > 1
            G = eye(3*Nstates)
        else
            G = eye(2*Nstates)
        end
        
        for iterscf = 1:MaxInnerSCF
            #
            V_loc = Ham.potentials.Hartree + Ham.potentials.XC
            if Ham.pspotNL.NbetaNL > 1
                VY = op_V_Ps_nloc( Ham, Y ) + op_V_loc( pw, V_loc, Y )
            else
                VY = op_V_loc( Ham.pw, V_loc, Y )
            end
            #
            A = real( T + Y'*VY )
            A = 0.5*( A + A' )
            #
            D, G = eig( A, B )
            λ = D[1:Nstates]
            #
            #E_kin = trace( G[:,set1]' * T * G[:,set1] )
            #println("Ekin from trace = ", E_kin)
            psi = Y*G[:,set1]
            rhoe = calc_rhoe( pw, Focc, psi )
            
            update!( Ham, rhoe )
            
        end

        # Calculate energies
        Energies = calc_energies( Ham, psi )
        Ham.energies = Energies
        Etot = Energies.Total
        diffE = abs( Etot - Etot_old )

        @printf("DCM: %5d %18.10f %18.10e\n", iter, Etot, diffE)

        if diffE < 1e-7
            @printf("DCM is converged: iter: %d , diffE = %10.7e\n", iter, diffE)
            break
        end

        Etot_old = Etot

        # No need to update potential, it is already updated in inner SCF loop

        # step
        if iter > 1
            P = Y[:,set4]*G[set4,set1]
        else
            P = Y[:,set2]*G[set2,set1]
        end


    end  # end of DCM iteration
    
    return λ, psi

end


