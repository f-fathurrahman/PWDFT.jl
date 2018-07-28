function KS_solve_DCM!( Ham::Hamiltonian;
                        NiterMax = 100, startingwfc=nothing,
                        savewfc=false, ETOT_CONV_THR=1e-6 )


	pw = Ham.pw
    Ngw = pw.gvecw.Ngw
    Ns = pw.Ns
    Npoints = prod(Ns)
    ΔV = pw.Ω/Npoints
    Focc = Ham.electrons.Focc
    Nstates = Ham.electrons.Nstates
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nspin = Ham.electrons.Nspin

    Nkspin = Nkpt*Nspin

    psiks = Array{Array{ComplexF64,2},1}(undef,Nkspin)

    #
    # Initial wave function
    #
    if startingwfc == nothing
        srand(1234)
        for ispin = 1:Nspin
        for ik = 1:Nkpt
            ikspin = ik + (ispin - 1)*Nkpt
            psiks[ikspin] = ortho_gram_schmidt( rand(ComplexF64,Ngw[ik],Nstates) )
        end
        end
    else
        psiks = startingwfc
    end

    #
    # Calculated electron density from this wave function and update Hamiltonian
    #
    Rhoe = zeros(Float64,Npoints,Nspin)
    for ispin = 1:Nspin
        idxset = (Nkpt*(ispin-1)+1):(Nkpt*ispin)
        Rhoe[:,ispin] = calc_rhoe( pw, Focc[:,idxset], psiks[idxset] )
    end
    update!(Ham, Rhoe)

    Rhoe_old = copy(Rhoe)

    evals = zeros(Float64,Nstates,Nkspin)
    
    # Starting eigenvalues and psi
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        Ham.ik = ik
        Ham.ispin = ispin
        ikspin = ik + (ispin - 1)*Nkpt
        evals[:,ikspin], psiks[ikspin] =
        diag_LOBPCG( Ham, psiks[ikspin], verbose_last=false, NiterMax=10 )
    end
    end

    Ham.energies = calc_energies( Ham, psiks )
    
    Etot_old = sum(Ham.energies)

    # subspace
    Y = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    R = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    P = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    G = Array{Array{ComplexF64,2},1}(undef,Nkspin)
    T = Array{Array{Float64,2},1}(undef,Nkspin)
    B = Array{Array{Float64,2},1}(undef,Nkspin)
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        ikspin = ik + (ispin - 1)*Nkpt
        Y[ikspin] = zeros( ComplexF64, Ngw[ik], 3*Nstates )
        R[ikspin] = zeros( ComplexF64, Ngw[ik], Nstates )
        P[ikspin] = zeros( ComplexF64, Ngw[ik], Nstates )
        G[ikspin] = zeros( ComplexF64, 3*Nstates, 3*Nstates )
        T[ikspin] = zeros( Float64, 3*Nstates, 3*Nstates )
        B[ikspin] = zeros( Float64, 3*Nstates, 3*Nstates )
    end
    end

    set1 = 1:Nstates
    set2 = Nstates+1:2*Nstates
    set3 = 2*Nstates+1:3*Nstates
    set4 = Nstates+1:3*Nstates

    MaxInnerSCF = 3

    for iter = 1:NiterMax
        
        for ispin = 1:Nspin
        for ik = 1:Nkpt
            Ham.ik = ik
            Ham.ispin = ispin
            ikspin = ik + (ispin - 1)*Nkpt
            #
            Hpsi = op_H( Ham, psiks[ikspin] )
            #
            psiHpsi = psiks[ikspin]' * Hpsi
            psiHpsi = 0.5*( psiHpsi + psiHpsi' )
            # Calculate residual
            R[ikspin] = Hpsi - psiks[ikspin]*psiHpsi
            R[ikspin] = Kprec( ik, pw, R[ikspin] )
            # Construct subspace
            Y[ikspin][:,set1] = psiks[ikspin]
            Y[ikspin][:,set2] = R[ikspin]
            #
            if iter > 1
                Y[ikspin][:,set3] = P[ikspin]
            end
            
            #
            # Project kinetic and ionic potential
            #
            KY = op_K( Ham, Y[ikspin] ) + op_V_Ps_loc( Ham, Y[ikspin] )
            #
            T[ikspin] = real( Y[ikspin]'*KY )
            B[ikspin] = real( Y[ikspin]'*Y[ikspin] )
            B[ikspin] = 0.5*( B[ikspin] + B[ikspin]' )

            if iter > 1
                G[ikspin] = Matrix(1.0I, 3*Nstates, 3*Nstates) #eye(3*Nstates)
            else
                G[ikspin] = Matrix(1.0I, 2*Nstates, 2*Nstates)
            end
        end
        end
        
        @printf("DCM iter: %3d\n", iter)

        for iterscf = 1:MaxInnerSCF
            
            for ispin = 1:Nspin
            for ik = 1:Nkpt
                #
                Ham.ik = ik
                Ham.ispin = ispin
                ikspin = ik + (ispin - 1)*Nkpt
                #
                V_loc = Ham.potentials.Hartree + Ham.potentials.XC[:,ispin]
                #
                if Ham.pspotNL.NbetaNL > 0
                    VY = op_V_Ps_nloc( Ham, Y[ikspin] ) + op_V_loc( ik, pw, V_loc, Y[ikspin] )
                else
                    VY = op_V_loc( ik, pw, V_loc, Y[ikspin] )
                end
                #
                A = real( T[ikspin] + Y[ikspin]'VY )
                A = 0.5*( A + A' )
                #
                D, G[ikspin] = eigen( A, B[ikspin] )
                evals[:,ikspin] = D[1:Nstates]
                #
                # update
                psiks[ikspin] = Y[ikspin]*G[ikspin][:,set1]
            end
            end
            

            for ispin = 1:Nspin
                idxset = (Nkpt*(ispin-1)+1):(Nkpt*ispin)
                Rhoe[:,ispin] = calc_rhoe( pw, Focc[:,idxset], psiks[idxset] )
            end
            update!( Ham, Rhoe )

            Rhoe_old = copy(Rhoe)

            # Calculate energies once again
            Ham.energies = calc_energies( Ham, psiks )
            Etot = sum(Ham.energies)
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
        Ham.energies = calc_energies( Ham, psiks )
        Etot = sum(Ham.energies)
        diffE = abs( Etot - Etot_old )
        @printf("DCM: %5d %18.10f %18.10e\n", iter, Etot, diffE)

        if abs(diffE) < ETOT_CONV_THR
            @printf("DCM is converged: iter: %d , diffE = %10.7e\n", iter, diffE)
            break
        end

        Etot_old = Etot

        # No need to update potential, it is already updated in inner SCF loop
        for ispin = 1:Nspin
        for ik = 1:Nkpt
            ikspin = ik + (ispin - 1)*Nkpt
            if iter > 1
                P[ikspin] = Y[ikspin][:,set4]*G[ikspin][set4,set1]
            else
                P[ikspin] = Y[ikspin][:,set2]*G[ikspin][set2,set1]
            end
        end
        end
    end  # end of DCM iteration
    
    Ham.electrons.ebands = evals[:,:]

    if savewfc
        for ikspin = 1:Nkpt*Nspin
            wfc_file = open("WFC_ikspin_"*string(ikspin)*".data","w")
            write( wfc_file, psiks[ikspin] )
            close( wfc_file )
        end
    end

    return
end


