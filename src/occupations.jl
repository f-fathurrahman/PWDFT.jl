"""
    calc_Focc(Nelectrons, wk, kT, evals, Nspin)

Determine occupation numbers for a given set of Kohn-Sham eigenvalues.
"""
function calc_Focc(
    Nelectrons::Float64,
    wk::Array{Float64,1},
    kT::Float64,
    evals::Array{Float64,2},
    Nspin::Int64;
    noncollinear = false
)
    @assert Nspin in [1,2]
    E_fermi = find_E_fermi( Nelectrons, wk, kT, evals, Nspin, noncollinear = noncollinear )
  
    Nkspin = size(evals)[2]
    Nstates = size(evals)[1]
  
    Focc = zeros(Nstates,Nkspin)

    for ikspin = 1:Nkspin, ist = 1:Nstates
        Focc[ist,ikspin] = wgauss( (E_fermi - evals[ist,ikspin])/kT )
    end
  
    if Nspin == 1 && !noncollinear
        return 2.0*Focc, E_fermi
    else
        return Focc, E_fermi
    end
  
end

function sum_kFocc( wk, Focc, Nspin )
    Nstates = size(Focc)[1]
    Nkspin = size(Focc)[2]

    wks = repeat(wk, Nspin)

    s = 0.0
    for ikspin = 1:Nkspin, ist = 1:Nstates
        s = s + wks[ikspin]*Focc[ist,ikspin]
    end
    return s
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



"""
    wgauss(x)

Compute Fermi distribution.
"""
function wgauss( x::Float64 )
# x = E_f - evals
# only n=-99 is implemented

    MAXARG = 200.0
    if x < -MAXARG
        return 0.0
    elseif x > MAXARG
        return 1.0
    else
        return 1.0/(1.0 + exp(-x))
    end

end


"""
    sumkg(wk, kT, evals, ene, Nspin)

Sum occupation numbers for a given set of Kohn-Sham eigenvalues up to `ene`.
"""
function sumkg(
  wk::Array{Float64,1},
  kT::Float64,
  evals::Array{Float64,2},
  ene::Float64,
  Nspin::Int64,
  noncollinear::Bool
)

    Nkspin = size(evals)[2]
    Nstates = size(evals)[1]

    wks = repeat(wk,Nspin)
    if (Nspin == 1) && !noncollinear
        wks[:] = wks[:]*2
    end

    ss = 0.0
    for ikspin = 1:Nkspin
        ss1 = 0.0
        for ist = 1:Nstates
            ss1 = ss1 + wgauss( (ene - evals[ist,ikspin])/kT )
        end
        ss = ss + wks[ikspin]*ss1
    end

    return ss
end

"""
    find_E_fermi(Nelectrons, wk, kT, evals, Nspin, [NiterMax=300, verbose=false])

Find E_fermi for a given set of Kohn-Sham eigenvalues.
"""
function find_E_fermi(
  Nelectrons::Float64,
  wk::Array{Float64,1},
  kT::Float64,
  evals::Array{Float64,2},
  Nspin::Int64;
  NiterMax=300, verbose=false, noncollinear = false
)

    Nkspin = size(evals)[2]
    Nstates = size(evals)[1]

    # determine lower and upper bound for bisection
    Elw = evals[1,1]
    Eup = evals[Nstates,1]
    for ikspin = 1:Nkspin
        Elw = min( Elw, evals[1,ikspin] )
        Eup = max( Eup, evals[Nstates,ikspin] )
    end
    Elw = Elw - 2*kT
    Eup = Eup + 2*kT

    sumklw = sumkg( wk, kT, evals, Elw, Nspin, noncollinear )
    sumkup = sumkg( wk, kT, evals, Eup, Nspin, noncollinear )

    SMALL = 1e-10

    if ( (sumkup - Nelectrons) < -eps() ) ||
       ( (sumklw - Nelectrons) >  eps() )
        error("Bounds for E_fermi is not found")
    end

    Ef = 0.5*(Eup + Elw)
    Ef_old = Ef
    for iter = 1:NiterMax
        sumkmid = sumkg( wk, kT, evals, Ef, Nspin, noncollinear )
        if abs(sumkmid-Nelectrons) < SMALL
            return Ef
        elseif sumkmid-Nelectrons < -SMALL
            Elw = Ef
        else
            Eup = Ef
        end
        Ef = 0.5*(Eup + Elw)
        diff_Ef = abs(Ef-Ef_old)
        if verbose
            @printf("find_E_fermi: %3d %18.10f %18.10f %18.10e\n", iter, Ef, sumkmid, diff_Ef)
        end
        if diff_Ef < SMALL
            return Ef
        end
        Ef_old = Ef
    end

    @printf("WARNING: Ef is not found after %d iterations\n", NiterMax)
    return Ef
    
end


"""
    w1gauss(x)

Calculate integ( y*delta(y,n), [y,-∞,x] )
where delta(y,n) is the approximation for the delta function.

# Notes

Adapted from file Modules/w1gauss.f90 from Quantum ESPRESSO. Only partially
implemented (for n=-99 only).
"""
function w1gauss(x)
    # n = -99 other cases are not yet implemented
    if abs(x) <= 36.0
        f = 1.0/(1.0 + exp(-x))
        onemf = 1.0 - f
        return f*log(f) + onemf*log(onemf)
    else
        return 0.0
    end
end


"""
    calc_entropy(wk, kT, evals, E_fermi, Nspin)

Calculate -TS term of the Kohn-Sham total energy.

# Arguments

- `wk`: weight for kpoints
- `kT`: smearing parameter
- `evals`: eigenvalues
- `E_fermi`: Fermi energy
- `Nspin`: number of spin component(s)
"""
function calc_entropy(
  wk::Array{Float64,1},
  kT::Float64,
  evals::Array{Float64,2},
  E_fermi::Float64,
  Nspin::Int64;
  noncollinear = false
)    
    Nkspin = size(evals)[2]
    Nstates = size(evals)[1]

    wks = repeat(wk,Nspin)
    if Nspin == 1 && !noncollinear
        wks[:] = wks[:]*2
    end

    mTS = 0.0
    for ikspin = 1:Nkspin
        for ist = 1:Nstates
            mTS = mTS + wks[ikspin]*kT*w1gauss( (E_fermi - evals[ist,ikspin])/kT )
        end
    end
    return mTS
end

function sum_upto_E_fermi(
    Focc::Array{Float64,2},
    evals::Array{Float64,2},
    wk::Array{Float64,1},
    E_fermi::Float64,
    Nspin::Int64
)
    Nstates = size(Focc)[1]
    Nkspin = size(Focc)[2]
    Nkpt = round(Int64, Nkspin/Nspin)

    sFocc = 0.0
    for ispin = 1:Nspin
    for ik = 1:Nkpt
        ikspin = ik + (ispin-1)*Nspin
        ss = 0.0
        for ist = 1:Nstates
            if evals[ist,ikspin] <= E_fermi
                ss = ss + Focc[ist,ikspin]
            end
        end
        sFocc = sFocc + wk[ik]*ss
    end
    end
    return sFocc
end

function sum_upto_E_fermi( Focc::Array{Float64,1}, evals::Array{Float64,1}, E_fermi::Float64 )
    sFocc = 0.0
    Nstates = size(Focc)[1]
    for ist = 1:Nstates
        if evals[ist] <= E_fermi
            sFocc = sFocc + Focc[ist]
        end
    end
    return sFocc
end

