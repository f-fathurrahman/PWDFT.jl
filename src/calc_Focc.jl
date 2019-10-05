function calc_Focc(
    Nelectrons::Float64,
    wk::Array{Float64,1},
    kT::Float64,
    evals::Array{Float64,2},
    Nspin::Int64
)

    E_fermi = find_E_fermi( Nelectrons, wk, kT, evals, Nspin )
  
    Nkspin = size(evals)[2]
    Nstates = size(evals)[1]
  
    Focc = zeros(Nstates,Nkspin)

    for ikspin = 1:Nkspin, ist = 1:Nstates
        Focc[ist,ikspin] = wgauss( (E_fermi - evals[ist,ikspin])/kT )
    end
  
    if Nspin == 1
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
