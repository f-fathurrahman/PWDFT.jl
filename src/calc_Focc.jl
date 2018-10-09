function calc_Focc( evals::Array{Float64,2}, wk::Array{Float64,1},
                    Nelectrons::Float64, kT::Float64;
                    Nspin=1, verbose=false )

    Nstates = size(evals)[1]
    Nkspin = size(evals)[2]
    Nkpt = round(Int64,Nkspin/Nspin)

    TOL = 1e-10
    MAXITER = 500

    Focc = zeros(Nstates,Nkspin)
    Nocc = max( round(Int64,Nelectrons/2), 1 ) # guard against Nelectrons==1
    if verbose
        @printf("Nelectrons = %d\n", Nelectrons)
        @printf("Nocc = %d\n", Nocc)
    end

    Focc_lb = zeros(Nstates,Nkspin)
    Focc_ub = zeros(Nstates,Nkspin)

    # use bisection to find E_fermi such that 
    #  sum_i Focc(i) = Nelectrons
    if Nstates > Nocc
        ilb = max(Nocc - 1,1)
        iub = Nocc + 1
        
        # FIXME: Need to guard against the case Nocc == 1
        # FIXME: spin-dependent ?
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
            flb = flb + sum( Focc_lb[:,ikspin] )*wk[ik]
            fub = fub + sum( Focc_ub[:,ikspin] )*wk[ik]
        end # ik
        end # ispin

        while ( (Nelectrons-flb)*(fub-Nelectrons) < 0 )
            if verbose
                @printf("calc_Focc: initial bounds are off:\n");
                @printf("flb = %18.10f, fub = %18.10f, Nelectrons = %d\n", flb, fub, Nelectrons)
            end
            if (flb > Nelectrons)
                if (ilb >= 1)
                    ilb = max(ilb - 1,1)
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
                    error("calc_Focc: cannot find a lower bound for E_fermi")
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
                    error("calc_Focc: cannot find an upper bound for E_fermi\n
                           Try increasing number of bands.")
                end
            end
            if ilb == 1
                break
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
        
        E_fermi_old = E_fermi

        for iter = 1:MAXITER
            diffNelec = abs(occsum-Nelectrons)

            if verbose
                @printf("%3d %18.10f %18.10f %18.10e\n", iter, E_fermi, occsum, diffNelec)
            end

            if diffNelec < TOL
                if verbose
                    @printf("Found E_fermi = %18.10f, occsum = %18.10f\n", E_fermi, occsum )
                end
                return Focc, E_fermi
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

            if abs(E_fermi - E_fermi_old) < eps()
                @printf("WARNING: Diff E_fermi is very small, bisection stops\n")
                return Focc, E_fermi
            end

            E_fermi_old = E_fermi

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

function check_Focc( Focc::Array{Float64,1}; Nspin=1 )
    if Nspin == 1
        for focc in Focc
            if focc > 2.0
                error("focc larger that 2.0 is found")
            end
        end
    else
        for focc in Focc
            if focc > 1.0
                error("focc larger that 1.0 is found")
            end
        end
    end
end
