
# spin-unpolarized version
function calc_Focc( evals::Array{Float64,2}, wk::Array{Float64,1},
                    Nelectrons::Float64, kT::Float64;
                    is_spinpol=false, verbose=false )

    Nstates = size(evals)[1]
    Nkpt = size(evals)[2]

    const TOL = 1e-10
    const MAXITER = 100

    Focc = zeros(Nstates,Nkpt)
    Nocc = round(Int64,Nelectrons/2)  # normally
    if verbose
        @printf("Nelectrons = %d\n", Nelectrons)
        @printf("Nocc = %d\n", Nocc)
    end

    Focc_lb = zeros(Nstates,Nkpt)
    Focc_ub = zeros(Nstates,Nkpt)

    # use bisection to find E_fermi such that 
    #  sum_i Focc(i) = Nelectrons
    if Nstates > Nocc
        ilb = Nocc - 1
        iub = Nocc + 1
        # FIXME: Need to guard against the case Nocc == 1
        lb = minimum(evals[ilb,:])
        ub = maximum(evals[iub,:])
        # make sure flb < Nelectrons and fub > Nelectrons
        # Use lb and ub as guess interval for searching Fermi energy
        flb = 0.0
        fub = 0.0
        for ik = 1:Nkpt
            Focc_lb[:,ik] = smear_FD(evals[:,ik], lb, kT, is_spinpol=is_spinpol)
            Focc_ub[:,ik] = smear_FD(evals[:,ik], ub, kT, is_spinpol=is_spinpol)
            #
            flb = flb + sum(Focc_lb[:,ik])*wk[ik]
            fub = fub + sum(Focc_ub[:,ik])*wk[ik]
        end
        #
        while ( (Nelectrons-flb)*(fub-Nelectrons) < 0 )
            if verbose
                @printf("calc_Focc: initial bounds are off:\n");
                @printf("flb = %18.10f, fub = %18.10f, Nelectrons = %d\n", flb, fub, Nelectrons)
            end
            if (flb > Nelectrons)
                if (ilb > 1)
                    ilb = ilb - 1
                    lb = minimum(evals[ilb,:])
                    flb = 0.0
                    for ik = 1:Nkpt
                        Focc[:,ik] = smear_FD(evals[:,ik], lb, kT, is_spinpol=is_spinpol)
                        flb = flb + sum(Focc[:,ik])*wk[ik]
                    end
                else
                    @printf("ERROR in calc_Focc: cannot find a lower bound for E_fermi\n")
                    exit()
                end
            end
            #
            if (fub < Nelectrons)
                if (iub < Nstates)
                    iub = iub + 1
                    ub  = maximum(evals[iub,:])
                    fub = 0.0
                    for ik = 1:Nkpt
                        Focc[:,ik] = smear_FD(evals[:,ik], ub, kT, is_spinpol=is_spinpol)
                        fub = fub + sum(Focc[:,ik])*wk[ik]
                    end
                else
                    @printf("ERROR in calc_Focc: cannot find an upper bound for E_fermi\n")
                    @printf("Try increasing the number of states\n")
                    exit()
                end
            end
        end  # while
        
        if verbose
            @printf("\nInitial bounds are OK: flb = %18.10f, fub = %18.10f\n", flb, fub)
        end
        
        E_fermi = (lb + ub)/2
        occsum = 0.0
        for ik = 1:Nkpt
            Focc[:,ik] = smear_FD(evals[:,ik], E_fermi, kT, is_spinpol=is_spinpol)
            occsum = occsum + sum(Focc[:,ik])*wk[ik]
        end
        
        for iter = 1:MAXITER
            diffNelec = abs(occsum-Nelectrons)
            if diffNelec < TOL
                if verbose
                    @printf("Found E_fermi = %18.10f, occsum = %18.10f\n", E_fermi, occsum )
                end
                return Focc, E_fermi
            end
            if verbose
                @printf("%3d %18.10f %18.10f %18.10e\n", iter, E_fermi, occsum, diffNelec)
            end
            if (occsum < Nelectrons)
                lb = E_fermi
            else
                ub = E_fermi
            end
            E_fermi = (lb + ub)/2
            occsum = 0.0
            for ik = 1:Nkpt
                Focc[:,ik] = smear_FD(evals[:,ik], E_fermi, kT, is_spinpol=is_spinpol)
                occsum = occsum + sum(Focc[:,ik])*wk[ik]
            end

        end #
        @printf("WARNING: Maximum iteration achieved, E_fermi is not found within specified tolarance\n")
        return Focc, E_fermi
    
    # 
    elseif (Nstates == Nocc)
        @printf("calc_Focc: Nstates is equal to Nocc\n")
        if is_spinpol
            for ik = 1:Nkpt
                Focc[:,ik] = 2.0*ones(Nstates)
            end
            E_fermi = maximum(evals[Nstates,:])
        else
            for ik = 1:Nkpt
                Focc[:,ik] = ones(Nstates)
            end
            E_fermi = maximum(evals[Nstates,:])
        end
        return Focc, E_fermi
    
    else
        @printf("ERROR: The number of eigenvalues in evals should be larger than Nelectrons")
        exit()
    end

end



# spin-unpolarized version
function calc_Focc( evals::Array{Float64,1}, Nelectrons::Float64, kT::Float64;
                    is_spinpol=false, verbose=false )

    Nstates = length(evals)
    const TOL = 1e-10
    const MAXITER = 100

    Focc = zeros(Nstates)
    Nocc = round(Int64,Nelectrons/2)  # normally
    if verbose
        @printf("Nelectrons = %d\n", Nelectrons)
        @printf("Nocc = %d\n", Nocc)
    end

    # use bisection to find E_fermi such that 
    #  sum_i Focc(i) = Nelectrons
    if Nstates > Nocc
        ilb = Nocc - 1
        iub = Nocc + 1
        # FIXME: Need to guard against the case Nocc == 1
        lb = evals[ilb]
        ub = evals[iub]
        # make sure flb < Nelectrons and fub > Nelectrons
        # Use lb and ub as guess interval for searching Fermi energy
        Focc_lb = smear_FD(evals, lb, kT, is_spinpol=is_spinpol)
        Focc_ub = smear_FD(evals, ub, kT, is_spinpol=is_spinpol)
        flb = sum(Focc_lb)
        fub = sum(Focc_ub)
        while ( (Nelectrons-flb)*(fub-Nelectrons) < 0 )
            if verbose
                @printf("calc_Focc: initial bounds are off:\n");
                @printf("flb = %18.10f, fub = %18.10f, Nelectrons = %d\n", flb, fub, Nelectrons)
            end
            if (flb > Nelectrons)
                if (ilb > 1)
                    ilb = ilb - 1
                    lb = evals[ilb]
                    Focc = smear_FD(evals, lb, kT, is_spinpol=is_spinpol)
                    flb = sum(Focc)
                else
                    @printf("ERROR in calc_Focc: cannot find a lower bound for E_fermi\n")
                    exit()
                end
            end
            #
            if (fub < Nelectrons)
                if (iub < Nstates)
                    iub = iub + 1
                    ub  = evals[iub]
                    Focc = smear_FD(evals, ub, kT, is_spinpol=is_spinpol)
                    fub = sum(Focc)
                else
                    @printf("ERROR in calc_Focc: cannot find an upper bound for E_fermi\n")
                    @printf("Try increasing the number of states\n")
                    exit()
                end
            end
        end  # while
        
        if verbose
            @printf("\nInitial bounds are OK: flb = %18.10f, fub = %18.10f\n", flb, fub)
        end
        
        E_fermi = (lb + ub)/2
        Focc = smear_FD(evals, E_fermi, kT, is_spinpol=is_spinpol)
        occsum = sum(Focc)
        
        for iter = 1:MAXITER
            diffNelec = abs(occsum-Nelectrons)
            if diffNelec < TOL
                if verbose
                    @printf("Found E_fermi = %18.10f, occsum = %18.10f\n", E_fermi, occsum )
                end
                return Focc, E_fermi
            end
            if verbose
                @printf("%3d %18.10f %18.10f %18.10e\n", iter, E_fermi, occsum, diffNelec)
            end
            if (occsum < Nelectrons)
                lb = E_fermi
            else
                ub = E_fermi
            end
            E_fermi = (lb + ub)/2
            Focc = smear_FD(evals, E_fermi, kT, is_spinpol=is_spinpol)
            occsum = sum(Focc)
        end #
        @printf("WARNING: Maximum iteration achieved, E_fermi is not found within specified tolarance\n")
        return Focc, E_fermi
    
    # 
    elseif (Nstates == Nocc)
        @printf("calc_Focc: Nstates is equal to Nocc\n")
        if is_spinpol
            Focc = 2.0*ones(Nstates)
            E_fermi = evals[Nstates]
        else
            Focc    = ones(Nstates)
            E_fermi = evals[Nstates]
        end
        return Focc, E_fermi
    
    else
        @printf("ERROR: The number of eigenvalues in evals should be larger than Nelectrons")
        exit()
    end

end
