function test_ElecVars( Ham::Hamiltonian )

    kT = 0.01
    evars = ElecVars(Ham)
    println(evars)

    E_fermi, mTS = update_occ!( Ham, evars, kT )
    
    println("E_fermi = ", E_fermi)

    println(Ham.electrons)

    println(Ham.electrons.ebands[:,1])

    Nspin = Ham.electrons.Nspin
    Nkspin = length(evars.psiks)
    Nkpt = Ham.pw.gvecw.kpoints.Nkpt
    Nstates = Ham.electrons.Nstates

    fprime = zeros(Float64,Nstates)
    fprimeNum = zeros(Float64,Nstates)
    dmuNum = zeros(Float64,Nspin)
    dmuDen = zeros(Float64,Nspin)
    
    #dmuNum[sIndex] += w * trace(fprime * ( diag(eVars.Hsub[q]) - eVars.Haux_eigs[q]) );
    #dmuDen[sIndex] += w * trace(fprime);

    w = copy(Ham.pw.gvecw.kpoints.wk)
    if Nspin == 1
        w = w*2.0
    end

    for ispin in 1:Nspin, ik in 1:Nkpt
        i = ik + (ispin - 1)*Nkpt
        ff = 0.0
        for ist in 1:Nstates
            fprime[ist] = smear_fermi_prime( evars.Haux_eigs[ist,i], E_fermi, kT )
            fprimeNum[ist] = fprime[ist] * ( real(evars.Hsub[i][ist,ist]) - evars.Haux_eigs[ist,i] )
        end
        #println(fprime)
        #println(fprimeNum)
        dmuNum[ispin] = dmuNum[ispin] + w[ik] * sum(fprimeNum)
        dmuDen[ispin] = dmuDen[ispin] + w[ik] * sum(fprime)
    end

    println("dmuNum = ", dmuNum)
    println("dmuDen = ", dmuDen)

    dmuContrib = sum(dmuNum)/sum(dmuDen)
    dBzContrib = 0.0

    println("dmuContrib = ", dmuContrib)

    # gradient w.r.t fillings except for constraint contributions
    # matrix gradF0 = eVars.Hsub[q] - eVars.Haux_eigs[q];
    gradF0 = zeros(ComplexF64,Nstates,Nstates)
    gradF = zeros(ComplexF64,Nstates,Nstates)
    #// gradient w.r.t fillings
    #matrix gradF = gradF0 - eye(eInfo.nBands)*eInfo.muEff(dmuContrib,dBzContrib,q);
      
    #grad->Haux[q] = qnum.weight * dagger_symmetrize(eInfo.smearGrad(eInfo.muEff(mu,Bz,q), eVars.Haux_eigs[q], gradF));
    
    #inline double muEff(double mu, double Bz, int q) const { return mu + Bz*qnums[q].spin; }
    #if( Kgrad ) //Drop the fermiPrime factors in preconditioned gradient:
    #  Kgrad->Haux[q] = (-e.cntrl.subspaceRotationFactor) * gradF0;


    g_Haux = deepcopy(evars.Hsub)
    Kg_Haux = deepcopy(evars.Hsub)

    g_tmp = zeros(ComplexF64,Nstates,Nstates)
    for ispin in 1:Nspin, ik in 1:Nkpt
        #
        i = ik + (ispin - 1)*Nkpt
        #
        gradF0[:] = evars.Hsub[i] - diagm( 0 => evars.Haux_eigs[:,i] )
        gradF[:] = copy(gradF0)
        for ist in 1:Nstates
            gradF[ist,ist] = gradF0[ist,ist] - dmuContrib # FIXME: not tested for spinpol
        end
        g_tmp[:] = grad_smear( smear_fermi, smear_fermi_prime, evars.Haux_eigs[:,i], E_fermi, kT, gradF )
        g_Haux[i] = w[ik] * 0.5 * (g_tmp' + g_tmp)
        #g_Haux[i] = 0.5 * (g_tmp' + g_tmp)
        #g_Haux[i] = copy(g_tmp)
        #g_Haux[i] = copy(gradF)
        Kg_Haux[i][:] = -copy(gradF0)
    end

    #
    println("g_Haux real")
    display(real(g_Haux[1])); println()
    println("g_Haux imag")
    display(imag(g_Haux[1])); println()
    #
    println("Kg_Haux real")
    display(real(Kg_Haux[1])); println()
    println("Kg_Haux imag")
    display(imag(Kg_Haux[1])); println()

    println("Pass here in test_ElecVars")
end