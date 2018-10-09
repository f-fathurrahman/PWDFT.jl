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

