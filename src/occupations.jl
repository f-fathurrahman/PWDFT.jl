
# x = E_f - evals
function wgauss( x::Float64; n=-99 )

    MAXARG = 200.0

    if n == -99
        if x < -MAXARG
            return 0.0
        elseif x > MAXARG
            return 1.0
        else
            return 1.0/(1.0 + exp(-x))
        end
    end

end


function sumkg(
  wk::Array{Float64,1},
  kT::Float64,
  evals::Array{Float64,2},
  ene::Float64,
  Nspin::Int64
)

    Nkspin = size(evals)[2]
    Nstates = size(evals)[1]

    wks = repeat(wk,Nspin)
    if Nspin == 1
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

function efermig(
  Nelectrons::Float64,
  wk::Array{Float64,1},
  kT::Float64,
  evals::Array{Float64,2},
  Nspin::Int64;
  NiterMax=300, verbose=false
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

    sumklw = sumkg( wk, kT, evals, Elw, Nspin )
    sumkup = sumkg( wk, kT, evals, Eup, Nspin )

    SMALL = 1e-10

    if ( (sumkup - Nelectrons) < -eps() ) ||
       ( (sumklw - Nelectrons) >  eps() )
        error("Bounds for E_fermi is not found")
    end

    Ef = 0.5*(Eup + Elw)
    Ef_old = Ef
    for iter = 1:NiterMax
        sumkmid = sumkg( wk, kT, evals, Ef, Nspin )
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
            @printf("efermig: %3d %18.10f %18.10f %18.10e\n", iter, Ef, sumkmid, diff_Ef)
        end
        if diff_Ef < SMALL
            return Ef
        end
        Ef_old = Ef
    end

    @printf("WARNING: Ef is not found after %d iterations\n", NiterMax)
    return Ef
    
end


function calc_Focc_v2(
  wk::Array{Float64,1},
  kT::Float64,
  evals::Array{Float64,2},
  E_fermi::Float64,
  Nspin::Int64
)

    Nkspin = size(evals)[2]
    Nstates = size(evals)[1]

    Focc = zeros(Nstates,Nkspin)

    for ikspin = 1:Nkspin
        for ist = 1:Nstates
            Focc[ist,ikspin] = wgauss( (E_fermi - evals[ist,ikspin])/kT )
        end
    end

    if Nspin == 1
        return 2.0*Focc
    else
        return Focc
    end

end
