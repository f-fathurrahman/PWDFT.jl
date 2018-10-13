"""
w1gauss(x,n) = integ( y*delta(y,n), [y,-∞,x] )
where delta(y,n) is the current approximation for the
delta function as obtained from w0gauss(x,n)

Adapted from file Modules/w1gauss.f90 from Quantum ESPRESSO.
"""
function w1gauss(x; n=-99)
    if n == -99
        if abs(x) <= 36.0
            f = 1.0/(1.0 + exp(-x))
            onemf = 1.0 - f
            return f*log(f) + onemf*log(onemf)
        else
            return 0.0
        end
    end
end


function calc_entropy_v2( wk, kT, evals, E_fermi, Nspin )
    
    Nkspin = size(evals)[2]
    Nstates = size(evals)[1]

    wks = repeat(wk,Nspin)
    if Nspin == 1
        wks[:] = wks[:]*2
    end

    mTS = 0.0
    for ikspin = 1:Nkspin
        for ist = 1:Nstates
            mTS = mTS + wks[ikspin]*kT*w1gauss( (E_fermi - evals[ist,ikspin])/kT )
            #println("wks = ", wks[ikspin])
            #println("mTS = ", mTS)
        end
    end
    return mTS
end
