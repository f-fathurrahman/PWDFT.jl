function my_dcm!( Ham )

    psiks = rand_BlochWavefunc(Ham)
    
    Rhoe = calc_rhoe(Ham, psiks)
    update!(Ham, Rhoe)
    
    Rhoe_new = similar(Rhoe)
    
    Ham.energies.NN = calc_E_NN(Ham.atoms)
    Ham.energies.PspCore = calc_PspCore_ene(Ham.atoms, Ham.pspots)

    Etot_old = 0.0
    Nconv = 0

    @assert length(psiks) == 1

    Nstates = Ham.electrons.Nstates

    set1 = 1:Nstates
    set2 = Nstates+1:2*Nstates
    set3 = 2*Nstates+1:3*Nstates
    set4 = Nstates+1:3*Nstates
    set5 = 1:2*Nstates

    psi = psiks[1]
    Nbasis = size(psi, 1)
    println("Nbasis = ", Nbasis)
    
    Hpsi = zeros(ComplexF64, Nbasis, Nstates)
    R = zeros(ComplexF64, Nbasis, Nstates)

    Y = zeros(ComplexF64, Nbasis, 3*Nstates)
    KY = zeros(ComplexF64, Nbasis, 3*Nstates)
    B = zeros(Float64, 3*Nstates, 3*Nstates)

    T = zeros(ComplexF64, 3*Nstates, 3*Nstates )

    ispin = 1
    V_loc = Ham.potentials.Hartree + Ham.potentials.XC[:,ispin]

    for iterDCM in 1:1

        Hpsi = Ham*psi
        Θ = psi' * Hpsi
        R = Kprec( 1, Ham.pw, Hpsi - psi*Θ )

        # Construct subspace
        Y[:,set1] = psi # XXX Probably use views?
        Y[:,set2] = R
        if iterDCM > 1
            Y[:,set3] = P[i]
        end

        if iterDCM > 1
            B = real(Y'*Y)
        else
            B[set5,set5] = real(Y[:,set5]' * Y[:,set5])
        end
        B = 0.5*(B + B')
        display(B)
        println()

        # Project Kinetic, Ps loc, and Ps nloc
        if iterDCM > 1
            KY = op_K( Ham, Y ) + op_V_Ps_loc( Ham, Y )
            if Ham.pspotNL.NbetaNL > 0
                KY = KY + op_V_Ps_nloc( Ham, Y )
            end
            T = real( Y' * KY )
        else
            # only set5=1:2*Nstates is active for iter=1
            KY[:,set5] = op_K( Ham, Y[:,set5] ) + op_V_Ps_loc( Ham, Y[:,set5] )
            if Ham.pspotNL.NbetaNL > 0
                KY[:,set5] = KY[:,set5] + op_V_Ps_nloc( Ham, Y[:,set5] )
            end
            T[set5,set5] = real( Y[:,set5]' * KY[:,set5] )
        end

        Rhoe_old = copy(Ham.rhoe)

        Etot_innerscf = 0.0
        Etot_innerscf_old = 0.0
        dEtot_innerscf = 1.0

    end

    println("my_dcm is finished")

    return

end

function DCM_inner_SCF!( iterDCM::Int64, Ham, KY, B, G, MaxInnerSCF::Int64 )

    for iterscf = 1:MaxInnerSCF
        
        for ispin = 1:Nspin, ik = 1:Nkpt
            #
            Ham.ik = ik
            Ham.ispin = ispin
            i = ik + (ispin - 1)*Nkpt
            #
            if iterDCM > 1
                yy = Y[i]
            else
                yy = Y[i][:,set5]
            end
            #
            if Ham.pspotNL.NbetaNL > 0
                VY = op_V_Ps_nloc( Ham, yy ) + op_V_loc( ik, pw, V_loc, yy )
            else
                VY = op_V_loc( ik, pw, V_loc, yy )
            end
            #
            if iter > 1
                A[i] = real( T[i] + yy'*VY )
                A[i] = 0.5*( A[i] + A[i]' )
            else
                aa = real( T[i][set5,set5] + yy'*VY )
                A[i] = 0.5*( aa + aa' )
            end
            #
            if iter > 1
                BG = B[i]*G[i][:,1:Nocc]
                C[i] = real( BG*BG' )
                C[i] = 0.5*( C[i] + C[i]' )
            else
                BG = B[i][set5,set5]*G[i][set5,1:Nocc]
                cc = real( BG*BG' )
                C[i][set5,set5] = 0.5*( cc + cc' )
            end
            #
            if iter > 1
                D[:,i], G[i] = eigen( A[i], B[i] )
            else
                D[set5,i], G[i][set5,set5] =
                eigen( A[i][set5,set5], B[i][set5,set5] )
            end
        #
        # update wavefunction
        if iter > 1
            psiks[i] = Y[i]*G[i][:,set1]
            ortho_sqrt!(psiks[i])  # is this necessary ?
        else
            psiks[i] = Y[i][:,set5]*G[i][set5,set1]
            ortho_sqrt!(psiks[i])
        end
        end # ispin, ik

    end


end