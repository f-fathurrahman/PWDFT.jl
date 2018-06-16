function KS_solve_DCM!( Ham::PWHamiltonian;
                        NiterMax = 100, startingwfc=nothing,
                        savewfc=true, ETOT_CONV_THR=1e-6 )


	pw = Ham.pw
    Ngw = pw.gvecw.Ngw
    Ns = pw.Ns
    Npoints = prod(Ns)
    ΔV = pw.Ω/Npoints
    Focc = Ham.electrons.Focc
    Nstates = Ham.electrons.Nstates
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt

    psik = Array{Array{ComplexF64,2},1}(undef,Nkpt)

    #
    # Initial wave function
    #
    if startingwfc == nothing
        srand(1234)
        srand(1234)
        for ik = 1:Nkpt
            psi = rand(Ngw[ik],Nstates) + im*rand(Ngw[ik],Nstates)
            psik[ik] = ortho_gram_schmidt(psi)
        end
    else
        psik = startingwfc
    end

    #
    # Calculated electron density from this wave function and update Hamiltonian
    #
    rhoe = calc_rhoe( pw, Focc, psik )
    update!(Ham, rhoe)

    rhoe_old = copy(rhoe)

    evals = zeros(Float64,Nstates,Nkpt)
    # Starting eigenvalues and psi
    for ik = 1:Nkpt
        Ham.ik = ik
        evals[:,ik], psik[ik] = diag_lobpcg( Ham, psik[ik], verbose_last=false, maxit=10 )
    end

    energies = calc_energies( Ham, psik )
    Ham.energies = energies
    Etot_old = energies.Total

    # subspace
    Y = Array{Array{ComplexF64,2},1}(undef,Nkpt)
    R = Array{Array{ComplexF64,2},1}(undef,Nkpt)
    P = Array{Array{ComplexF64,2},1}(undef,Nkpt)
    G = Array{Array{ComplexF64,2},1}(undef,Nkpt)
    T = Array{Array{Float64,2},1}(undef,Nkpt)
    B = Array{Array{Float64,2},1}(undef,Nkpt)
    for ik = 1:Nkpt
        Y[ik] = zeros( ComplexF64, Ngw[ik], 3*Nstates )
        R[ik] = zeros( ComplexF64, Ngw[ik], Nstates )
        P[ik] = zeros( ComplexF64, Ngw[ik], Nstates )
        G[ik] = zeros( ComplexF64, 3*Nstates, 3*Nstates )
        T[ik] = zeros( Float64, 3*Nstates, 3*Nstates )
        B[ik] = zeros( Float64, 3*Nstates, 3*Nstates )
    end

    set1 = 1:Nstates
    set2 = Nstates+1:2*Nstates
    set3 = 2*Nstates+1:3*Nstates
    set4 = Nstates+1:3*Nstates

    MaxInnerSCF = 3

    for iter = 1:NiterMax
        
        for ik = 1:Nkpt
            Ham.ik = ik
            Hpsi = op_H( Ham, psik[ik] )
            #
            psiHpsi = psik[ik]' * Hpsi
            psiHpsi = 0.5*( psiHpsi + psiHpsi' )
            #
            R[ik] = Hpsi - psik[ik]*psiHpsi
            R[ik] = Kprec( ik, pw, R[ik] )
            #
            Y[ik][:,set1] = psik[ik]
            Y[ik][:,set2] = R[ik]
            if iter > 1
                Y[ik][:,set3] = P[ik]
            end
            
            #
            # Project kinetic and ionic potential
            #
            KY = op_K( Ham, Y[ik] ) + op_V_Ps_loc( Ham, Y[ik] )
            #
            T[ik] = real( Y[ik]'*KY )
            B[ik] = real( Y[ik]'*Y[ik] )
            B[ik] = 0.5*( B[ik] + B[ik]' )

            if iter > 1
                G[ik] = Matrix(1.0I, 3*Nstates, 3*Nstates) #eye(3*Nstates)
            else
                G[ik] = Matrix(1.0I, 2*Nstates, 2*Nstates)
            end
        end
        
        @printf("DCM iter: %3d\n", iter)

        for iterscf = 1:MaxInnerSCF
            ispin = 1  # FIXME setup properly for loop over Nspin
            V_loc = Ham.potentials.Hartree + Ham.potentials.XC[:,ispin]

            for ik = 1:Nkpt
                Ham.ik = ik
                if Ham.pspotNL.NbetaNL > 0
                    VY = op_V_Ps_nloc( Ham, Y[ik] ) + op_V_loc( ik, pw, V_loc, Y[ik] )
                else
                    VY = op_V_loc( ik, pw, V_loc, Y[ik] )
                end
                #
                A = real( T[ik] + Y[ik]'*VY )
                A = 0.5*( A + A' )
                #
                D, G[ik] = eigen( A, B[ik] )
                evals[:,ik] = D[1:Nstates]
                #
                
                #E_kin = trace( G[:,set1]' * T * G[:,set1] )
                #println("Ekin from trace = ", E_kin)

                # update psik
                psik[ik] = Y[ik]*G[ik][:,set1]
            end
            
            #rhoe = 0.9*calc_rhoe( pw, Focc, psi ) + 0.1*rhoe_old # use mixing
            rhoe = calc_rhoe( pw, Focc, psik )
            
            update!( Ham, rhoe )

            rhoe_old = copy(rhoe)

            # Calculate energies once again
            Ham.energies = calc_energies( Ham, psik )
            Etot = Ham.energies.Total
            diffE = -(Etot - Etot_old)

            @printf("innerSCF: %5d %18.10f %18.10e", iterscf, Etot, diffE)
            # positive value of diffE is taken as reducing
            if diffE < 0.0
                @printf(" : Energy is not reducing !\n")
            else
                @printf("\n")
            end
            
        end

        # Calculate energies once again
        energies = calc_energies( Ham, psik )
        Ham.energies = energies
        Etot = energies.Total
        diffE = abs( Etot - Etot_old )
        @printf("DCM: %5d %18.10f %18.10e\n", iter, Etot, diffE)

        if abs(diffE) < ETOT_CONV_THR
            @printf("DCM is converged: iter: %d , diffE = %10.7e\n", iter, diffE)
            break
        end

        Etot_old = Etot

        # No need to update potential, it is already updated in inner SCF loop

        for ik = 1:Nkpt
            if iter > 1
                P[ik] = Y[ik][:,set4]*G[ik][set4,set1]
            else
                P[ik] = Y[ik][:,set2]*G[ik][set2,set1]
            end
        end
    end  # end of DCM iteration
    
    Ham.electrons.ebands = evals[:,:]

    if savewfc
        for ik = 1:Nkpt
            wfc_file = open("WFC_k_"*string(ik)*".data","w")
            write( wfc_file, psik[ik] )
            close( wfc_file )
        end
    end

    return
end


