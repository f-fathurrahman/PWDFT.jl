
# spin-unpolarized version
function calc_Focc( evals::Array{Float64,2}, wk::Array{Float64,1},
                    Nelectrons::Float64, kT::Float64;
                    Nspin=1, verbose=false )

    Nstates = size(evals)[1]
    Nkspin = size(evals)[2]
    Nkpt = round(Int64,Nkspin/Nspin)
    
    println("Nstates = ", Nstates)
    println("Nkspin = ", Nkspin)
    println("Nkpt = ", Nkpt)

    const TOL = 1e-10
    const MAXITER = 100

    Focc = zeros(Nstates,Nkspin)
    Nocc = round(Int64,Nelectrons/2)  # normally
    if verbose
        @printf("Nelectrons = %d\n", Nelectrons)
        @printf("Nocc = %d\n", Nocc)
    end

    Focc_lb = zeros(Nstates,Nkspin)
    Focc_ub = zeros(Nstates,Nkspin)

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
        for ispin = 1:Nspin
        for ik = 1:Nkpt
            ikspin = ik + (ispin - 1)*Nkpt
            Focc_lb[:,ikspin] = smear_FD(evals[:,ikspin], lb, kT, Nspin=Nspin)
            Focc_ub[:,ikspin] = smear_FD(evals[:,ikspin], ub, kT, Nspin=Nspin)
            #
            flb = flb + sum(Focc_lb[:,ikspin])*wk[ik]
            fub = fub + sum(Focc_ub[:,ikspin])*wk[ik]
        end # ik
        end # ispin
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
                    for ispin = 1:Nspin
                    for ik = 1:Nkpt
                        ikspin = ik + (ispin - 1)*Nkpt
                        Focc[:,ikspin] = smear_FD(evals[:,ikspin], lb, kT, Nspin=Nspin)
                        flb = flb + sum(Focc[:,ikspin])*wk[ik]
                    end # ik
                    end # ispin
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
                    for ispin = 1:Nspin
                    for ik = 1:Nkpt
                        ikspin = ik + (ispin - 1)*Nkpt
                        Focc[:,ikspin] = smear_FD(evals[:,ikspin], ub, kT, Nspin=Nspin)
                        fub = fub + sum(Focc[:,ikspin])*wk[ik]
                    end # ik
                    end # ispin
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
        for ispin = 1:Nspin
        for ik = 1:Nkpt
            ikspin = ik + (ispin - 1)*Nkpt
            Focc[:,ikspin] = smear_FD(evals[:,ikspin], E_fermi, kT, Nspin=Nspin)
            occsum = occsum + sum(Focc[:,ikspin])*wk[ik]
        end
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
            for ispin = 1:Nspin
            for ik = 1:Nkpt
                ikspin = ik + (ispin - 1)*Nkpt
                Focc[:,ikspin] = smear_FD(evals[:,ikspin], E_fermi, kT, Nspin=Nspin)
                occsum = occsum + sum(Focc[:,ikspin])*wk[ik]
            end
            end

        end #
        @printf("WARNING: Maximum iteration achieved, E_fermi is not found within specified tolarance\n")
        return Focc, E_fermi
    
    # 
    elseif (Nstates == Nocc)
        @printf("calc_Focc: Nstates is equal to Nocc\n")
        #
        if Nspin == 2
            # spin-polarized
            for ispin = 1:Nspin
            for ik = 1:Nkpt
                ikspin = ik + (ispin - 1)*Nkpt
                Focc[:,ikspin] = ones(Nstates)
            end
            end
            E_fermi = maximum(evals[Nstates,:])
        else
            # No spin-polarization
            for ispin = 1:Nspin
            for ik = 1:Nkpt
                ikspin = ik + (ispin - 1)*Nkpt
                Focc[:,ik] = 2.0*ones(Nstates)
            end
            end
            E_fermi = maximum(evals[Nstates,:])
        end
        return Focc, E_fermi
    
    else
        @printf("ERROR: The number of eigenvalues in evals should be larger than Nelectrons")
        exit()
    end

end



# spin-unpolarized version, wk is 1
# UNUSED ??
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
