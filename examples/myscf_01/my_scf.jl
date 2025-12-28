function my_scf!(
    Ham::Hamiltonian, psiks::BlochWavefunc;
    startingwfc=:random,
    startingrhoe=:gaussian,
    betamix=0.2, NiterMax=100, verbose=true,
    check_rhoe=false,
    use_smearing=false, kT=1e-3,
    update_psi="LOBPCG",
    etot_conv_thr=1e-6, cheby_degree=8
)

    pw = Ham.pw
    Ngw = pw.gvecw.Ngw
    wk = Ham.pw.gvecw.kpoints.wk
    #
    kpoints = pw.gvecw.kpoints
    Nkpt = kpoints.Nkpt
    #
    Ns = pw.Ns
    Npoints = prod(Ns)
    CellVolume = pw.CellVolume
    dVol = CellVolume/Npoints
    #
    electrons = Ham.electrons
    Nelectrons = electrons.Nelectrons
    Focc = copy(electrons.Focc) # make sure to use the copy
    Nstates = electrons.Nstates
    Nstates_occ = electrons.Nstates_occ
    Nspin = electrons.Nspin_channel
    Nkspin = Nkpt*Nspin
    Nstates_occ = electrons.Nstates_occ

    #
    # Calculated electron density from this wave function and update Hamiltonian
    #
    Rhoe = zeros(Float64,Npoints,Nspin)

    if startingrhoe == :gaussian && startingwfc == :random
        Rhoe = guess_rhoe_atomic( Ham )
    else
        calc_rhoe!(Ham, psiks, Rhoe)
    end

    update!(Ham, psiks, Rhoe)

    Etot_old = 0.0
    Rhoe_new = zeros(Float64,Npoints,Nspin)
    diffRhoe = zeros(Nspin)
    evals = zeros(Float64,Nstates,Nkspin)
    ETHR_EVALS_LAST = 1e-8

    ethr = 0.1

    @printf("\n")
    @printf("Self-consistent iteration begins ...\n")
    @printf("update_psi = %s\n", update_psi)
    @printf("Density mixing with betamix = %10.5f\n", betamix)
    if use_smearing
        @printf("Smearing = %f\n", kT)
    end
    println("") # blank line before SCF iteration info

    # calculate E_NN
    Ham.energies.NN = calc_E_NN( Ham.atoms )

    CONVERGED = 0

    Etot = 0.0
    Etot_old = Inf

    @assert Nspin == 1  # Current limitation

    # Preallocated arrays
    epsxc = zeros(Float64, Npoints)
    Vin = zeros(Float64, Npoints)

    for iterSCF in 1:NiterMax

        #@printf("ethr = %18.10e\n", ethr)
        evals, ethr = do_diag!(Ham, psiks, update_psi, iterSCF,
            Nstates_occ, ethr, ETHR_EVALS_LAST, cheby_degree)

        if use_smearing
            Focc, E_fermi = calc_Focc( Nelectrons, wk, kT, evals, Nspin )
            Ham.energies.mTS = calc_entropy( wk, kT, evals, E_fermi, Nspin )
            Ham.electrons.Focc = copy(Focc)
        end

        calc_rhoe!(Ham, psiks, Rhoe_new)
        for ispin in 1:Nspin
            @views diffRhoe[ispin] = norm(Rhoe_new[:,ispin] - Rhoe[:,ispin])
        end

        Eband = 0.0
        for ispin in 1:Nspin, ik = 1:Nkpt
            i = ik + (Nspin - 1)*Nkpt
            for ist in 1:Nstates
                Eband = Eband + Focc[ist,i]*wk[ik]*evals[ist,i]
            end
        end
        
        @views Vin[:] = Ham.potentials.Hartree[:] + Ham.potentials.XC[:]
        @views deband = -dot( Vin, Rhoe_new[:,1] )*dVol

        # Mix
        @views Rhoe[:] = betamix*Rhoe_new[:] + (1-betamix)*Rhoe[:]

        update!( Ham, psiks, Rhoe )

        descf = 0.0
        for ip in 1:Npoints
           dρ = Rhoe[ip,1] - Rhoe_new[ip,1]
           V = Ham.potentials.Hartree[ip] + Ham.potentials.XC[ip,1]
           descf += V*dρ
        end
        descf *= dVol

        @views EHartree = 0.5*dot(Ham.potentials.Hartree, Rhoe[:,1])*dVol

        if Ham.xcfunc == "PBE"
            calc_epsxc_PBE!( Ham.xc_calc, Ham.pw, Ham.rhoe, epsxc )
        else
            if Ham.rhoe_core == nothing
                calc_epsxc_VWN!( Ham.xc_calc, Ham.rhoe, epsxc )
            else
                calc_epsxc_VWN!( Ham.xc_calc, Ham.rhoe + Ham.rhoe_core )
            end
        end
        if Ham.rhoe_core == nothing
            Exc = dot( epsxc, Ham.rhoe ) * dVol
        else
            Exc = dot( epsxc, Ham.rhoe + Ham.rhoe_core ) * dVol
        end

        if use_smearing

        end
        
        #Etot = sum(Ham.energies)
        Etot = Eband + deband + EHartree + Exc + Ham.energies.NN + descf + Ham.energies.mTS
        diffE = abs( Etot - Etot_old )

        @printf("SCF: %8d %18.10f %18.10e %18.10e\n", iterSCF, Etot, diffE, diffRhoe[1] )
        
        if diffE < etot_conv_thr
            CONVERGED = CONVERGED + 1
        else  # reset CONVERGED
            CONVERGED = 0
        end

        if CONVERGED >= 2
            @printf("SCF is converged: iterSCF: %d , diffE = %10.7e\n", iterSCF, diffE)
            println()
            println("Total energy terms in Ry:")
            @printf("Eband    = %18.10f Ry\n", 2*Eband)
            @printf("deband   = %18.10f Ry\n", 2*deband)
            @printf("OneEle   = %18.10f Ry\n", 2*(Eband + deband))
            @printf("EHartree = %18.10f Ry\n", 2*EHartree)
            @printf("Exc      = %18.10f Ry\n", 2*Exc)
            @printf("NN       = %18.10f Ry\n", 2*Ham.energies.NN)
            @printf("mTS      = %18.10f Ry\n", 2*Ham.energies.mTS)
            @printf("-----------------------------\n")
            @printf("Total    = %18.10f Ry\n", 2*Etot)
            break
        end
        #
        Etot_old = Etot

        flush(stdout)
    end

    Ham.electrons.ebands[:] = evals[:] # FIXME

    return
end



function do_diag!(
    Ham,
    psiks,
    update_psi,
    iter,
    Nstates_occ,
    ethr, ETHR_EVALS_LAST,
    cheby_degree
)

    # determine convergence criteria for diagonalization
    if iter == 1
        ethr = 0.1
    elseif iter == 2
        ethr = 0.01
    else
        ethr = ethr/5.0
        ethr = max( ethr, ETHR_EVALS_LAST )
    end

    if update_psi == "LOBPCG"
        evals =
        diag_LOBPCG!( Ham, psiks, verbose=false, verbose_last=false,
                      Nstates_conv=Nstates_occ )

    elseif update_psi == "davidson"
        evals =
        diag_davidson!( Ham, psiks, verbose=false, verbose_last=false,
                        Nstates_conv=Nstates_occ )                

    elseif update_psi == "PCG"
        evals =
        diag_Emin_PCG!( Ham, psiks, verbose=false, verbose_last=false,
                        Nstates_conv=Nstates_occ )

    elseif update_psi == "CheFSI"
        evals =
        diag_CheFSI!( Ham, psiks, cheby_degree )

    else
        error( @sprintf("Unknown method for update_psi = %s\n", update_psi) )
    end

    return evals, ethr

end

